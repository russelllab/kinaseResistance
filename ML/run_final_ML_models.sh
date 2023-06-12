# script to generate and save the best models and
# scalers for each predictor

# AIvLD
python ML.py --s AIvLD --m AIvLD --c columns_to_consider.txt 5 7 5 100 AIvLD

# AIvNLD
python ML.py --s AIvNLD --m AIvNLD --c columns_to_consider.txt 5 5 7 100 AIvNLD

# LDvNAI
python ML.py --s LDvNAI --m LDvNAI --c columns_to_consider.txt 5 10 7 100 LDvNAI

# RvN
python ML.py --s RvN --m RvN --c columns_to_consider.txt 5 5 3 100 RvN