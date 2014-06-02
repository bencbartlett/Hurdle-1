#Data Compiler
#Takes output from 31 different files and stiches them together.
compiledString = ""
for i in xrange(1,32):
    filename = "Final_Sphere_Heights_Composite_"+str(i)+".dat"
    filehandle = open(filename, "r")
    print "Reading ",filename
    string = filehandle.read()
    compiledString += string
    filehandle.close()
compiledString = "{"+compiledString.replace(" ,","")
print compiledString
filehandle = open("FINAL_END.dat", "w")
print "Writing data..."
filehandle.write(compiledString)
filehandle.close()