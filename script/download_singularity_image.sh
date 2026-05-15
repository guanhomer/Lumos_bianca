SINGULARITY_CACHE="/home/xkwash/Lumos/singularity"  
mkdir -p $SINGULARITY_CACHE  
  
# Pull all containers used by Lumos  
singularity pull \
  $SINGULARITY_CACHE/quay.io-jmonlong-minimap2_samtools-v2.24_v1.16.1.img \
  docker://quay.io/jmonlong/minimap2_samtools:v2.24_v1.16.1  

singularity pull \
  $SINGULARITY_CACHE/hkubal-clair3-v1.0.11.img \
  docker://hkubal/clair3:v1.0.11

singularity pull \
  $SINGULARITY_CACHE/mkolmogo-longphase-1.7.3.img \
  docker://mkolmogo/longphase:1.7.3

singularity pull \
  $SINGULARITY_CACHE/mkolmogo-whatshap-2.3.img \
  docker://mkolmogo/whatshap:2.3

singularity pull \
  $SINGULARITY_CACHE/mkolmogo-modkit-0.4.1.img \
  docker://mkolmogo/modkit:0.4.1

singularity pull \
  $SINGULARITY_CACHE/gokcekeskus-severus-v1_6.img \
  docker://gokcekeskus/severus:v1_6

singularity pull \
  $SINGULARITY_CACHE/mkolmogo-wakhan-0.4.0.img \
  docker://mkolmogo/wakhan:0.4.0

singularity pull \
  $SINGULARITY_CACHE/google-deepsomatic-1.9.0.img \
  docker://google/deepsomatic:1.9.0