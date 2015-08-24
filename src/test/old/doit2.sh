MAT=/acrm/usr/local/share/hb/data/hbmatrices.dat
PDB=/acrm/data/pdb/pdb1tsr.ent
CUT=0.2
./checkhbond -c $CUT -m $MAT $PDB B126 B131 ASN 
./checkhbond -c $CUT -m $MAT $PDB B131 B126 ASP 
./checkhbond -c $CUT -m $MAT $PDB B127 B282 ARG 
./checkhbond -c $CUT -m $MAT $PDB B127 B286 GLU 
./checkhbond -c $CUT -m $MAT $PDB B132 B271 GLU 
./checkhbond -c $CUT -m $MAT $PDB B132 B285 GLU 
./checkhbond -c $CUT -m $MAT $PDB B140 B198 GLU 
./checkhbond -c $CUT -m $MAT $PDB B146 B144 GLN 
./checkhbond -c $CUT -m $MAT $PDB B155 B259 ASP 
./checkhbond -c $CUT -m $MAT $PDB B158 B215 SER 
./checkhbond -c $CUT -m $MAT $PDB B158 B258 GLU 
./checkhbond -c $CUT -m $MAT $PDB B163 B249 ARG 
./checkhbond -c $CUT -m $MAT $PDB B249 B168 HIS 
./checkhbond -c $CUT -m $MAT $PDB B183 B175 ARG 
./checkhbond -c $CUT -m $MAT $PDB B183 B196 ARG 
./checkhbond -c $CUT -m $MAT $PDB B192 B207 ASP 
./checkhbond -c $CUT -m $MAT $PDB B214 B207 ASP 
./checkhbond -c $CUT -m $MAT $PDB B236 B253 THR 
./checkhbond -c $CUT -m $MAT $PDB B280 B281 ASP 
