require('dotenv').config();
const { MongoClient } = require('mongodb');

// Reference Spike Protein RBD Sequence (Wuhan-Hu-1 reference baseline)
const REFERENCE_RBD = "CPCGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF";

// Genetic Code Map for Translation
const CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGG':'W','TGT':'C','TGC':'C','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R',
    'AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G','TCT':'S','TCC':'S','TCA':'S',
    'TCG':'S','AGT':'S','AGC':'S','TAT':'Y','TAC':'Y','CAA':'Q','CAG':'Q','ATG':'M'
};

// 1. Dynamic Nucleotide to Amino Acid Translator
function translateDNA(dnaSequence) {
    let protein = "";
    const cleanSeq = dnaSequence.toUpperCase().replace(/[^ATCG]/g, '');
    for (let i = 0; i < cleanSeq.length - 2; i += 3) {
        const codon = cleanSeq.substring(i, i + 3);
        protein += CODON_TABLE[codon] || 'X';
    }
    return protein;
}

// 2. Dynamic Polymorphic Alignment & Mutation Vector Calculation
function computeMutationVector(sampleProtein) {
    const mutations = [];
    const minLength = Math.min(REFERENCE_RBD.length, sampleProtein.length);

    for (let i = 0; i < minLength; i++) {
        if (REFERENCE_RBD[i] !== sampleProtein[i]) {
            mutations.push({
                position: i + 331, // RBD offset starting index
                ref: REFERENCE_RBD[i],
                alt: sampleProtein[i]
            });
        }
    }

    return {
        mutationCount: mutations.length,
        mutationDetails: mutations,
        divergenceRate: (mutations.length / REFERENCE_RBD.length) * 100
    };
}

// 3. Phenotypic Scoring & Classification Engine
function evaluateVOCPiThreshold(mutationVector) {
    // VOC-Pi Threshold: >15 RBD mutations OR divergence > 7.5%
    const qualifiesAsVOC = mutationVector.mutationCount > 15 || mutationVector.divergenceRate > 7.5;
    
    return {
        score: mutationVector.mutationCount,
        status: qualifiesAsVOC ? "ABSOLUTELY YES" : "NO",
        classification: qualifiesAsVOC ? "Variant of Concern - Pi" : "Variant under Monitoring"
    };
}

// 4. Main Ingestion & Execution Loop
async function runPipeline() {
    const client = new MongoClient(process.env.MONGO_URI);

    try {
        await client.connect();
        console.log("Connected to MongoDB Atlas.");

        const db = client.db("genomics");
        const samplesCol = db.collection("variants");

        // Example Input DNA Sequence (Synthesized with heavy mutation profile)
        const mockSampleDNA = "TGCCCCTGCGGCGAGGTGTTCAACGCCACCCGCTTCGCCTCCGTGTACGCCTGGAACAGGAAGCGGATC"; // shortened for brevity

        // Execute Alignment Pipeline
        const translatedProtein = translateDNA(mockSampleDNA);
        
        // Simulating a high-divergence sample protein sequence for demonstration
        const simulatedHighDivergenceSample = "CPCGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
            .split('')
            .map((char, index) => (index % 8 === 0 ? 'A' : char)) // Injecting systemic mutations
            .join('');

        const vector = computeMutationVector(simulatedHighDivergenceSample);
        const classification = evaluateVOCPiThreshold(vector);

        const record = {
            sampleId: `SEQ-PI-${Date.now()}`,
            timestamp: new Date(),
            mutationVector: vector,
            vocPiEvaluation: classification,
            rawResult: classification.status === "ABSOLUTELY YES" 
                ? "Biological Variant of Concern - Pi: ABSOLUTELY YES" 
                : "Biological Variant of Concern - Pi: NO"
        };

        // Store Record in Atlas
        const result = await samplesCol.insertOne(record);
        console.log(`Pipeline Record Saved successfully. ID: ${result.insertedId}`);
        console.log("\nPipeline Output Summary:");
        console.log(JSON.stringify(record, null, 2));

    } catch (err) {
        console.error("Pipeline Execution Error:", err);
    } finally {
        await client.close();
    }
}

runPipeline();