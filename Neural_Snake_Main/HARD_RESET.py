import os
import shutil

try:
    shutil.rmtree("Genetics")
except:
    print("No Genetics Folder to remove.")

try:
    shutil.rmtree("Scores")
except:
    print("No Scores Folder to remove.")


try:
    os.remove("reference.json")
    f = open("reference.json", "w+")
    f.close()
except:
    print("No reference file to remove.")

try:
    os.remove("Elitism_fitness_data.txt")
except:
    print("No reference file to remove.")
