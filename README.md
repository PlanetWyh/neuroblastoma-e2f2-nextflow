# Neuroblastoma E2F2 Regulatory Network and Survival Analysis



## Overview



This project investigates the role of E2F2 transcriptional activity in neuroblastoma using both single-cell and bulk RNA sequencing datasets.



The project consists of two complementary analyses:



### Stage 1: Single-cell Neuroblastoma Analysis



Using the GSE216155 single-cell RNA-seq dataset, neuroendocrine (NE) tumor cells were identified and analyzed using pySCENIC to reconstruct transcription factor regulatory networks and quantify E2F2 regulon activity.



### Stage 2: Bulk Neuroblastoma Survival Analysis



Using TARGET neuroblastoma RNA-seq data, E2F2 target gene activity was quantified using ssGSEA. Associations between E2F2 activity and patient survival were evaluated using Kaplan–Meier survival analysis and log-rank testing. Findings were validated using the independent GSE49711 neuroblastoma cohort.



---



## Repository Structure



```text

neuroblastoma_nextflow/

├── main.nf

├── nextflow.config

├── scripts/

│   ├── stage1/

│   └── stage2/

│   └── requirements.txt

├── resources/

│   ├── pyscenic_resources/

├── data/

│   ├── stage1/

│   └── stage2/

├── results/

└── README.md

```



---



## Workflow



### Stage 1: Single-cell Analysis



```text

GSE216155

    ↓

load_nb.py

    ↓

preprocess_nb.py

    ↓

subset_ne.py

    ↓

clean_ne.py

    ↓

export_for_pyscenic.py

    ↓

GRNBoost2

    ↓

pySCENIC Context Analysis

    ↓

AUCell

    ↓

E2F2 Regulon Visualization

```



### Stage 2: Bulk RNA-seq Survival Analysis



```text

TARGET RNA-seq

    ↓

build_target_matrix.py

    ↓

calculate_e2f2_ssgsea.py

    ↓

prepare_survival_data.py

    ↓

survival_analysis.py

    ↓

validate_e2f2_gse49711.py

    ↓

plot_gse49711_e2f2.py

```



---



## Software Requirements



* Python 3.10+

* Nextflow 26+

* Java 11+



Python dependencies are listed in:



```text

requirements.txt

```



---



## Running the Pipeline



### Stage 1 only



```bash

nextflow run main.nf \

    --run_stage1 true \

    --run_stage2 false

```



### Stage 2 only



```bash

nextflow run main.nf \

    --run_stage1 false \

    --run_stage2 true

```



### Full Analysis



```bash

nextflow run main.nf \

    --run_stage1 true \

    --run_stage2 true

```



---



## Outputs



### Stage 1



* Processed neuroblastoma AnnData objects

* pySCENIC regulons

* AUCell regulon activity scores

* E2F2 regulon UMAP visualization



### Stage 2



* TARGET expression matrix

* E2F2 ssGSEA scores

* Survival analysis tables

* Kaplan–Meier plots

* Validation analyses using GSE49711



---



## Reproducibility



The complete workflow is implemented using Nextflow DSL2 to ensure reproducibility, automated execution, dependency tracking, and transparent generation of intermediate and final outputs.



## Data Availability



Due to file size limitations, raw datasets are not included in this repository. Data and resources can be downloaded [here](https://drive.google.com/drive/folders/1Ur0Cbm0ajBEcnn7C02Q97-qfm8TCOe31?usp=share_link).



Datasets used:



- GSE216155 (single-cell neuroblastoma)

- TARGET Neuroblastoma RNA-seq

- GSE49711 neuroblastoma cohort



Required files should be downloaded and placed into the `data/` directory according to the structure described in `nextflow.config`.
