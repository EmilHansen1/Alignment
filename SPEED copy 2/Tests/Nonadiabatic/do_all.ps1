Write-Host "Running Na2 Singlet calculation..."
python .\Na2Singlet.py

Write-Host "Running K2 Singlet calculations..."
python .\K2SingletHomo.py
python .\K2SingletHetero.py

Write-Host "Running Rb2 Singlet calculations..."
python .\Rb2SingletHomo.py
python .\Rb2SingletHetero.py

Write-Host "Running K2 Triplet calculations..."
python .\K2TripletHomo.py
python .\K2TripletHetero.py

Write-Host "Running Rb2 Triplet calculations..."
python .\Rb2TripletHomo.py
python .\Rb2TripletHetero.py