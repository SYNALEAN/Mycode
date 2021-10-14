import json
import sys
import gc

for chro in range(0,455):
    print('strat '+str(chro))
    with open(str(chro)+".txt", "r") as f:  # 打开文件
        data = f.readlines()
        for line in data:
            line = line.replace("\n", '').upper()  # 去掉列表中每一个元素的换行符 #去掉列表中每一个元素的换行符
        print(str(len(data)) + 'lines; finished reading')

    


    slicedlength = 11
    reference = {}

    for index in range(0, len(data[0]) - slicedlength):
        theKey = data[0][index:index + slicedlength]
        countn = theKey.find('a')+theKey.find('t')+theKey.find('c')+theKey.find('g')+theKey.find('N')
        #print(index)
        if countn==-5:
            #print(sys.getsizeof(reference))
            if theKey in reference:
                reference[theKey].append(index)
            else:
                reference[theKey] = [index]

    
    
    
    json_str = json.dumps(reference, indent=4)
    with open('./json/'+str(chro)+'_reference.json', 'w') as json_file:
        json_file.write(json_str)

    del data
    del reference
    gc.collect()
    print('clear')
