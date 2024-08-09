# config_and_print.py
import scipy.io

data_path = '../../projects/GSE63525/GM12878/'
resolutions = [1000000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
data_types = ['oe']
genomeID = "hg19"
hic_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"
dark_regions_hg19_url = 'https://www.encodeproject.org/files/ENCFF000EHJ/@@download/ENCFF000EHJ.bigWig'
mappability_threshold = 0.6
H3K4me3_GM12878_hg19_url = "https://www.encodeproject.org/files/ENCFF154XCY/@@download/ENCFF154XCY.bigWig"
iterations = 400

config = {
    'data_path': data_path,
    'resolutions': resolutions,
    'chromosomes': chromosomes,
    'data_types': data_types,
    'genomeID': genomeID,
    'hic_url': hic_url,
    'dark_regions_hg19_url': dark_regions_hg19_url,
    'mappability_threshold': mappability_threshold,
    'H3K4me3_GM12878_hg19_url': H3K4me3_GM12878_hg19_url,
    'iterations': iterations
}

scipy.io.savemat('config.mat', config)

print(f"data_path='{data_path}'")
print(f"resolutions=({', '.join(map(str, resolutions))})")
print(f"chromosomes=({' '.join(chromosomes)})")
print(f"data_types=({' '.join(data_types)})")
print(f"genomeID='{genomeID}'")
print(f"hic_url='{hic_url}'")
print(f"dark_regions_hg19_url='{dark_regions_hg19_url}'")
print(f"mappability_threshold={mappability_threshold}")
print(f"H3K4me3_GM12878_hg19_url={H3K4me3_GM12878_hg19_url}")
print(f"iterations={iterations}")
