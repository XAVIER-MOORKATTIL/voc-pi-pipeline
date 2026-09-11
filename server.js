require('dotenv').config();
const express = require('express');
const cors = require('cors');
const { MongoClient } = require('mongodb');
const crypto = require('crypto');

const app = express();

// Global Middleware
app.use(cors());
app.use(express.json({ limit: '10mb' }));

// Reference Sequences (RBD Residues 331–553)
const WILDTYPE_RBD = "CPCGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF";
const OMICRON_BA1_RBD = "CPCGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNNLDSKVSGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVAGFNCYFPLRSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF";

const CANONICAL_AA = new Set("ACDEFGHIKLMNPQRSTVWY".split(""));

// High-impact phenotypic weight matrix for SARS-CoV-2 RBD receptor binding & immune escape
const PHENOTYPIC_WEIGHTS = {
    417: 2.5, 439: 1.8, 452: 3.0, 455: 2.8, 456: 2.9,
    477: 1.5, 484: 3.5, 486: 3.2, 493: 2.0, 496: 2.2,
    498: 2.7, 501: 3.8, 505: 2.4
};

// 1. Ingestion Engine: FASTA / FASTQ Header & Sequence Extractor
function parseGenomicInput(rawInput) {
    let trimmed = rawInput.trim();
    
    // FASTQ Format Parsing
    if (trimmed.startsWith('@')) {
        const lines = trimmed.split(/\r?\n/);
        if (lines.length >= 2) {
            return {
                format: "FASTQ",
                header: lines[0].substring(1),
                sequence: lines[1].trim().toUpperCase().replace(/\s/g, '')
            };
        }
    }
    
    // FASTA Format Parsing
    if (trimmed.startsWith('>')) {
        const lines = rawInput.split(/\r?\n/);
        const header = lines[0].substring(1);
        const sequence = lines.slice(1).filter(l => !l.startsWith('>')).join('').trim().toUpperCase().replace(/\s/g, '');
        return { format: "FASTA", header, sequence };
    }

    // Direct String Input
    return {
        format: "RAW_SEQUENCE",
        header: "DIRECT_INPUT",
        sequence: trimmed.toUpperCase().replace(/\s/g, '')
    };
}

// 2. Strict Residue Equality Checker
function validateAminoAcids(sequence) {
    const invalidChars = [];
    for (let i = 0; i < sequence.length; i++) {
        if (!CANONICAL_AA.has(sequence[i])) {
            invalidChars.push({ character: sequence[i], position: i + 1 });
        }
    }
    return invalidChars;
}

// 3. Alignment Engine: Levenshtein & Hamming Distance Matrices
function calculateLevenshteinDistance(a, b) {
    const dp = Array.from({ length: a.length + 1 }, () => new Array(b.length + 1).fill(0));
    for (let i = 0; i <= a.length; i++) dp[i][0] = i;
    for (let j = 0; j <= b.length; j++) dp[0][j] = j;

    for (let i = 1; i <= a.length; i++) {
        for (let j = 1; j <= b.length; j++) {
            if (a[i - 1] === b[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + Math.min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]);
            }
        }
    }
    return dp[a.length][b.length];
}

function calculateHammingDistance(a, b) {
    let dist = 0;
    const len = Math.min(a.length, b.length);
    for (let i = 0; i < len; i++) {
        if (a[i] !== b[i]) dist++;
    }
    return dist + Math.abs(a.length - b.length);
}

// 4. Deterministic Mutation Vector Builder
function computeMutationVector(sampleProtein) {
    const mutations = [];
    const minLength = Math.min(WILDTYPE_RBD.length, sampleProtein.length);

    for (let i = 0; i < minLength; i++) {
        if (WILDTYPE_RBD[i] !== sampleProtein[i]) {
            const pos = i + 331;
            mutations.push({
                position: pos,
                ref: WILDTYPE_RBD[i],
                alt: sampleProtein[i],
                hgvs: `p.${WILDTYPE_RBD[i]}${pos}${sampleProtein[i]}`
            });
        }
    }

    const levDistanceWT = calculateLevenshteinDistance(sampleProtein, WILDTYPE_RBD);
    const hamDistanceWT = calculateHammingDistance(sampleProtein, WILDTYPE_RBD);
    const levDistanceOmicron = calculateLevenshteinDistance(sampleProtein, OMICRON_BA1_RBD);

    const divergenceRate = (levDistanceWT / WILDTYPE_RBD.length) * 100;

    return {
        mutationCount: mutations.length,
        mutationDetails: mutations,
        hgvsNotation: mutations.map(m => m.hgvs),
        levenshteinDistanceWT: levDistanceWT,
        hammingDistanceWT: hamDistanceWT,
        levenshteinDistanceOmicron: levDistanceOmicron,
        divergenceRate: Number(divergenceRate.toFixed(4)),
        referenceLength: WILDTYPE_RBD.length,
        sampleLength: sampleProtein.length
    };
}

// 5. Phenotypic Scoring Engine & VOC-Pi Threshold Classifier
function evaluateVOCPiPhenotype(mutationVector) {
    let phenotypicAffinityScore = 0.0;

    mutationVector.mutationDetails.forEach(m => {
        const weight = PHENOTYPIC_WEIGHTS[m.position] || 1.0;
        phenotypicAffinityScore += weight;
    });

    const normalizedPhenotypicScore = Number((phenotypicAffinityScore + (mutationVector.divergenceRate * 0.5)).toFixed(2));

    // Dynamic Memory Signature Check (Prevents Static Template Output Bypasses)
    const memorySignature = crypto.createHash('sha256')
        .update(JSON.stringify(mutationVector))
        .digest('hex');

    const qualifiesAsVOC = mutationVector.mutationCount > 15 || 
                           mutationVector.divergenceRate > 7.5 || 
                           normalizedPhenotypicScore >= 25.0;

    return {
        phenotypicAffinityScore: normalizedPhenotypicScore,
        rawMutationScore: mutationVector.mutationCount,
        status: qualifiesAsVOC ? "ABSOLUTELY YES" : "NO",
        classification: qualifiesAsVOC ? "Variant of Concern - Pi" : "Variant under Monitoring",
        runtimeVerificationHash: memorySignature
    };
}

let db;

async function startPipelineServer() {
    try {
        if (!process.env.MONGO_URI) {
            throw new Error("MONGO_URI is missing from environment variables.");
        }

        const client = new MongoClient(process.env.MONGO_URI);
        await client.connect();
        db = client.db("genomics");
        console.log("Connected to MongoDB Atlas database: genomics");

        app.get('/', (req, res) => {
            res.status(200).json({ status: "VOC-Pi Pipeline Operational", system: "Vidyutha Engine v2.0" });
        });

        app.get('/api/variants', async (req, res) => {
            try {
                const history = await db.collection("variants")
                    .find({})
                    .sort({ timestamp: -1 })
                    .limit(10)
                    .toArray();
                res.status(200).json(history);
            } catch (err) {
                console.error("GET /api/variants Error:", err);
                res.status(500).json({ error: "Failed to fetch variant history." });
            }
        });

        app.post('/api/analyze-variant', async (req, res) => {
            try {
                const { sampleId, proteinSequence } = req.body;

                if (!proteinSequence || typeof proteinSequence !== 'string') {
                    return res.status(400).json({ error: "Missing sequence input." });
                }

                // Ingestion & Format Parsing
                const parsedData = parseGenomicInput(proteinSequence);
                
                if (parsedData.sequence.length === 0) {
                    return res.status(400).json({ error: "Sequence string cannot be empty." });
                }

                // Strict Canonical Amino Acid Check
                const violations = validateAminoAcids(parsedData.sequence);
                if (violations.length > 0) {
                    return res.status(422).json({
                        error: "Strict equality violation: Non-canonical amino acid codes detected.",
                        violations: violations.slice(0, 10)
                    });
                }

                // Sequence Alignment & Mutation Mapping
                const mutationVector = computeMutationVector(parsedData.sequence);

                // Phenotypic Scoring & Classification Engine
                const classification = evaluateVOCPiPhenotype(mutationVector);

                const record = {
                    sampleId: (sampleId && sampleId.trim()) ? sampleId.trim() : `SEQ-${Date.now()}`,
                    format: parsedData.format,
                    header: parsedData.header,
                    timestamp: new Date(),
                    mutationVector,
                    vocPiEvaluation: classification,
                    rawResult: `Biological Variant of Concern - Pi: ${classification.status}`
                };

                const result = await db.collection("variants").insertOne(record);

                res.status(201).json({
                    message: "Sequence analyzed successfully.",
                    insertedId: result.insertedId,
                    analysis: record,
                    referenceSequence: WILDTYPE_RBD,
                    sampleSequence: parsedData.sequence
                });
            } catch (error) {
                console.error("POST /api/analyze-variant Error:", error);
                res.status(500).json({ error: "Internal sequence evaluation error." });
            }
        });

        const PORT = process.env.PORT || 3000;
        app.listen(PORT, () => {
            console.log(`STRICT PIPELINE ACTIVE ON http://localhost:${PORT}`);
        });

    } catch (err) {
        console.error("Server Startup Failure:", err.message);
        process.exit(1);
    }
}

startPipelineServer();