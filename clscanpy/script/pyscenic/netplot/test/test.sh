module load OESingleCell/3.0.d

# Rscript /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/oepySCENIC/oepySCENIC/script/netplot/main.r \
#   --input /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/test/DZOE2024091605/pySCENIC/data/expr_mat.adjacencies.tsv \
#   --interestTF E2f7,E2f1 \
#   --outdir /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/netplot/test \
#   --targetTop 100 \
#   --showTargetTop 0.3 \
#   --highlightTarget Syce2,Donson,Pdxp

Rscript /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/oepySCENIC/oepySCENIC/script/netplot/main.r \
  --input /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/test/DZOE2024091605/pySCENIC/data/expr_mat.adjacencies.tsv \
  --interestTF E2f7,E2f1 \
  --outdir /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/oepySCENIC/oepySCENIC/script/netplot/test

