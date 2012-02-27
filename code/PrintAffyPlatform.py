import glob,os,sys,time,posix
import utilities
from mycel import MyCEL

celFilePath = sys.argv[1]

cel = MyCEL()

print cel.read_platform(celFilePath)
