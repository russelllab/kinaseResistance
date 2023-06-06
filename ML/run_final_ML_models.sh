# script to generate and save the best models and
# scalers for each predictor

# AIvLD
python ML.py --s AIvLD --m AIvLD 5 3 5 50 AIvLD

# AIvNLD
python ML.py --s AIvNLD --m AIvNLD 5 4 4 100 AIvNLD

# LDvNAI
python ML.py --s LDvNAI --m LDvNAI 5 5 5 50 LDvNAI

# RvN
python ML.py --s RvN --m RvN 5 4 4 100 RvN