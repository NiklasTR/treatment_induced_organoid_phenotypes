#rsync -rav /Users/rindtorf/Documents/GitHub/promise/data/ rindtorf@b110-sc2cn01:/home/rindtorf/github/promise/data/ 
rsync -rav /rsession-store/rindtorf/promise rindtorf@b110-sc2cn01:/home/rindtorf/data/promise/github

module load R/3.6.2
Rscript /rsession-store/rindtorf/promise/make_vignettes.R

find ./* -size +100M | cat >> .gitignore
find ./* -size +100M | cat >> transfer_file