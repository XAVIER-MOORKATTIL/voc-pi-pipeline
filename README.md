# VOC-Pi Genomic Pipeline Dashboard

An automated SARS-CoV-2 genomic alignment and variant classification engine designed to detect post-Omicron "Variant of Concern - Pi" strains.

## Key Features

* **Flexible Sequence Ingestion**: Native parsing for raw sequence strings, FASTA (`>Header`), and FASTQ (`@Header`) formats.
* **Alignment Engine**: Levenshtein edit-distance and Hamming distance computations against Wildtype Spike RBD and Omicron BA.1 lineage profiles.
* **Deterministic Mutation Mapping**: Computes amino acid substitution vectors and transcribes exact HGVS notations (`p.X331Y`).
* **Phenotypic Affinity Scoring Engine**: Evaluates site-specific receptor binding domain (RBD) affinity weights across critical positions (e.g., 417, 452, 484, 501, 505).
* **Deterministic Anti-Bypass Security**: Real-time SHA-256 runtime memory signatures verify dynamic variable resolution to prevent static output bypasses.
* **Audit & Storage**: Direct integration with MongoDB Atlas for execution logging and CSV audit generation.

## Environment Variables

Create a `.env` file in the root directory:

```env
PORT=3000
MONGO_URI=mongodb+srv://<username>:<password>@cluster.mongodb.net/genomics?retryWrites=true&w=majority

Local Development Setup
Clone the repository:

Bash
git clone [https://github.com/your-username/voc-pi-pipeline.git](https://github.com/your-username/voc-pi-pipeline.git)
cd voc-pi-pipeline
Install dependencies:

Bash
npm install express cors mongodb dotenv
Run the server:

Bash
node server.js
Access the frontend dashboard:
Open index.html in your web browser or serve it via Live Server (http://127.0.0.1:5500/index.html).