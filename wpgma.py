#!usr/bin/python
import array
import numpy

#This function compares the two sequence strings and tell the which string is smaller and which is larger sequence

def calscore(str1, str2):
    itr = 0
    score = 0
    if len(str1) > len(str2):
        itr = len(str2)
        score = len(str1) - len(str2)
    else:
        itr = len(str1)
        score = len(str2) - len(str1)
    for i in range(itr):
        if str1[i] != str2[i]:
            score += 1
    return score
    
#This function take the input from the user and aligns in the matrix

def seqtoscore():
    filename = input("Enter the file name \n")
    seq = []
    with open(filename) as f:
        for ln in f:
            for x in ln.split(','):
                seq.append(x)


    score = [[0 for i in range(len(seq))] for j in range(len(seq))]

    for i in range(len(seq)):
        for j in range(len(seq)):
            if i == j :
                score[i][j]=0
            else:
                score[i][j]=calscore(seq[i],seq[j])


    return score
    


# This function reads the distance matrix from the file that user uploaded

def readInput():      
    matrix = []       
    matrix = seqtoscore()
    #print(matrix)
    matrix = numpy.asarray(matrix)
    matrix = matrix.astype(numpy.float)
    length = len(matrix)
    print ("Distance matrix:")
    print (matrix)
    print ("" )
    print ("Number of sequences:", length)
    print ("")
    return matrix, length

# This function finds the minimum value in the distance matrix

def matrixMinimum(matrix, length):
    min_index_i = 0
    min_index_j = 0
    minimum=float('inf')
    for i in range (length):

        temp = matrix[i]
        min_value = numpy.min(temp[numpy.nonzero(temp)])
        j = temp.tolist().index(min_value)

        if min_value < minimum: 
            min_index_i = i
            min_index_j = j
            minimum = min_value
   
    return min_index_i,min_index_j

def wpgma(matrix, length,dictionary):
    leaves = []
    count=0
    for i in range (0, length):    
        leaves.append("S"+str(i+1))
    
    numberOfClusters=i+1
    
    while(length>1):
        
        
        numberOfClusters=numberOfClusters+1
        count=count+1
        min_index_i,min_index_j = matrixMinimum(matrix, length)
        
        leaves.append("S"+str(numberOfClusters))  
        distance=matrix[min_index_i][min_index_j]/float(2)
        
        size=0
        if leaves[min_index_i] not in dictionary.keys():
            size1=1
            distance1=distance
        else:
            size1=dictionary[leaves[min_index_i]][4]
            distance1=distance-max(dictionary[leaves[min_index_i]][0],dictionary[leaves[min_index_i]][2])
        
        if leaves[min_index_j] not in dictionary.keys():
            size2=1
            distance2=distance
        else:
            size2=dictionary[leaves[min_index_j]][4]
            distance2=distance-max(dictionary[leaves[min_index_j]][0],dictionary[leaves[min_index_j]][2])
        dictionary["S"+str(numberOfClusters)]=[distance1,leaves[min_index_i],distance2,leaves[min_index_j],size1+size2]
        
        # Create a new row and column
        matrix = numpy.insert(matrix, length, values=float(0), axis=0)
        matrix = numpy.insert(matrix, length, values=float(0), axis=1)
        
        for i in range (0, length):
            matrix[-1][i]=matrix[i][-1] = (size1*matrix[i][min_index_i] + size2*matrix[i][min_index_j])/float(size1+size2)
        
        # Delete the minimum value
        if min_index_i < min_index_j:
            matrix = numpy.delete(matrix, min_index_i, 0)
            matrix = numpy.delete(matrix, min_index_i, 1)
            matrix = numpy.delete(matrix, (min_index_j)-1, 0)
            matrix = numpy.delete(matrix, (min_index_j)-1, 1)
            length = len(matrix)
            del leaves[min_index_j]
            del leaves[min_index_i]

        else:
            matrix = numpy.delete(matrix, min_index_i, 0)
            matrix = numpy.delete(matrix, min_index_i, 1)
            matrix = numpy.delete(matrix, min_index_j, 0)
            matrix = numpy.delete(matrix, min_index_j, 1)            
            length = len(matrix)
            del leaves[min_index_i]
            del leaves[min_index_j]
        
    return "S"+str(numberOfClusters)

def printCluster(dictionary,finalCluster):
    stack=[]
    result=[]
    stack.append(finalCluster)
    while stack:
        
        current=stack.pop()
        if isinstance( current , float ):
            if isinstance( current_prev , float ):
                result.pop()
                result.append(")")
            result.append(":"+str(current))
            result.append(",")
            
        elif current in dictionary.keys():
        
            stack.append(dictionary[current][0])
            stack.append(dictionary[current][1])
            stack.append(dictionary[current][2])
            stack.append(dictionary[current][3])
            result.append("(")
        else:
            result.append(current)
        current_prev=current
        
    result.pop()
    result.append(")")  
    return result

if __name__ == "__main__":
    matrix, length = readInput()
    dictionary={}
    finalCluster = wpgma(matrix, length,dictionary)
    result=printCluster(dictionary,finalCluster)
    result=''.join(result)
    result=result+";"

#ete3 is tool for pylogenetic tree construction
    
    from ete3 import Tree
    tree=Tree(result)
    print ("WPGMA Resultant Clustering:")
    print (""  )
    print (result)
    print (""  )
    print (tree)

