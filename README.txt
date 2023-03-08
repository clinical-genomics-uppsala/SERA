
      ___           ___           ___           ___     
     /\__\         /\__\         /\  \         /\  \    
    /:/ _/_       /:/ _/_       /::\  \       /::\  \   
   /:/ /\  \     /:/ /\__\     /:/\:\__\     /:/\:\  \  
  /:/ /::\  \   /:/ /:/ _/_   /:/ /:/  /    /:/ /::\  \ 
 /:/_/:/\:\__\ /:/_/:/ /\__\ /:/_/:/__/___ /:/_/:/\:\__\
 \:\/:/ /:/  / \:\/:/ /:/  / \:\/:::::/  / \:\/:/  \/__/
  \::/ /:/  /   \::/_/:/  /   \::/~~/~~~~   \::/__/     
   \/_/:/  /     \:\/:/  /     \:\~~\        \:\  \     
     /:/  /       \::/  /       \:\__\        \:\__\    
     \/__/         \/__/         \/__/         \/__/    


# Build singularity

```bash
# Build docker
docker build -f Dockerfile.sera .  --tag gmsuppsala/sera:1.0.0
docker save gmsuppsala/sera:1.0.0 -o gmsuppsala_sera_1.0.0.tar
singularity build gmsuppsala_sera_1.0.0.simg  docker-archive://gmsuppsala_sera_1.0.0.tar

docker build -f Dockerfile.ampliconMapping .  --tag gmsuppsala/ampliconmapping:0.0.1
docker save gmsuppsala/ampliconmapping:0.0.1 -o gmsuppsala_ampliconmapping_0.0.1.tar
singularity build gmsuppsala_ampliconmapping_0.0.1.simg  docker-archive://gmsuppsala_ampliconmapping_0.0.1.tar
```

# Start pipeline

## Create input data
```
singularity run -B /beegfs-storage -B /projects -B /data -B /opt /beegfs-storage/projects/wp4/nobackup/singularity_cache/gmsuppsala_sera_1.0.0.simg python3 /home/patsm159/SERA/bin/pythonscript/createInputFile_moriarty.py -a klinik -g MARVIN -i /projects/wp1/nobackup/ngs/klinik/sample_files/2031/20310310_PS_index.csv -n annovar -p wp1 -refDir /data/ref_data/wp1/refFiles_20230123/refFiles
```

## Run SERA
```
module load sera/VERSION
SERA_PBS_SUBMIT.sh -p /projects/wp1/nobackup/ngs/klinik/analys/2031/20310310_PS -i /projects/wp1/nobackup/ngs/klinik/analys/2031/20310310_PS/inputFile  -f
```