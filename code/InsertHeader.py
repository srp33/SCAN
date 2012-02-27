import os, sys, glob
import utilities

inFilePath = sys.argv[1]
header = sys.argv[2].decode('string-escape')

text = header + "\n" + utilities.readTextFromFile(inFilePath)
utilities.writeScalarToFile(text, inFilePath)
