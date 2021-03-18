
def readLogFile():
    mylines = []
    with open("ID001_001log.txt","rt") as myfile: 
        for myline in myfile:
            mylines.append(myline.rstrip('\n'))


    #Find receiver position
    substrR = "Receiver: "
    for line in mylines:          # string to be searched
        index = 0                   # current index: character being compared
        prev = 0                    # previous index: last character compared
        while index < len(line):    # While index has not exceeded string length,
            index = line.find(substrR, index)  # set index to first occurrence of substring
            if index == -1:           # If nothing was found,
                break                   # exit the while loop. 
            receiverPos = "(" + line[index+len(substrR):index+len(substrR)+17] + ")"
                                                 
            prev = index + len(substrR)       # remember this position for next loop.
            index += len(substrR)      # increment the index by the length of substr.
                                      # (Repeat until index > line length)
                              
    #Find source position
    substrS = "Source: "                  
    for line in mylines:          
        index = 0                 
        prev = 0                   
        while index < len(line):   
            index = line.find(substrS, index) 
            if index == -1:         
                break                 
            sourcePos = "(" + line[index+len(substrS):index+len(substrS)+17] + ")"
                                             
            prev = index + len(substrS)      
            index += len(substrS)      
                             
                              
    #Find water depth
    substrD = "Level: "                  
    for line in mylines:         
        index = 0                 
        prev = 0                    
        while index < len(line):   
            index = line.find(substrD, index) 
            if index == -1:          
                break                  
            waterDepth = line[index+len(substrD):index+len(substrD)+10]
                                             
            prev = index + len(substrD)     
            index += len(substrD)
                           
                              
    #Find sampling frequency
    substrF = "Frequency: "                  
    for line in mylines:         
        index = 0                 
        prev = 0                
        while index < len(line):   
            index = line.find(substrF, index) 
            endIndex = line.find("(Hz)")
            if index == -1:           
                break                   
            samplingFreq = line[index+len(substrF):endIndex] + "Hz"
                                             
            prev = index + len(substrF)      
            index += len(substrF)     
                             
                              
    #Find bandwidth
    substrBmin = "from "                  
    for line in mylines:          
        index = 0                 
        prev = 0                  
        while index < len(line):   
            index = line.find(substrBmin, index) 
            if index != 0:
                endIndex = line.find('.00')
            endIndex = line.find('Hz')
            if index == -1:           
                break                  
            bandwidth = line[index+len(substrBmin):endIndex] + " - " + line[index+len(substrBmin)+15:index+len(substrBmin)+28]
            fmin = line[index+len(substrBmin):index+len(substrBmin)+8]                  
               
            prev = index + len(substrBmin)       
            index += len(substrBmin)     
                          

    #Find water temperature
    substrT = "Temp: "                  
    for line in mylines:          
        index = 0                   
        prev = 0                  
        while index < len(line):    
            index = line.find(substrT, index) 
            if index == -1:         
                break                  
            waterTemp = line[index+len(substrT):index+len(substrT)+8]
                                             
            prev = index + len(substrT)       
            index += len(substrT)     
                             
                              
    #Find signal length
    substrL = "length: "                  
    for line in mylines:          
        index = 0                  
        prev = 0                    
        while index < len(line):   
            index = line.find(substrL, index) 
            endIndex = line.find("   Trailing")
            if index == -1:           # If nothing was found,
                break                   # exit the while loop. 
            signalLength = line[index+len(substrL):endIndex]
                                             
            prev = index + len(substrL)       # remember this position for next loop.
            index += len(substrL)      # increment the index by the length of substr.
                              # (Repeat until index > line length)
                              
    print(bandwidth,waterTemp,sourcePos,receiverPos,samplingFreq,waterDepth)
    return bandwidth,waterTemp,sourcePos,receiverPos,samplingFreq,signalLength,waterDepth,fmin
                              