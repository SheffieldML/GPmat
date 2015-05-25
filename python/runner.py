import reader as reader
from pylab import *
import time
from learning import nonLinearLearner as nLL
from learning import linearLearner as lL
from graphics import Graphics as G

class runner():
    def run(self,bvh):
        readInst = reader.MyReader(bvh,True)
        dt ,points, limits = readInst.read()

        #Create new linear learner object
        learner = lL(readInst.channels)

        #Create new non-linear learner object
        #nLLearner = nLL(readInst.channels)

        # Create a new graphics object (no animation)
        #G(dt,points,limits,learner,readInst)

        # Create a new animation (no learning)
        G(dt,points,limits,learner,readInst,True)

        #run the test method
        self.myTest(readInst.channels,nLLearner)

    # Split a set of channels into test data and training data
    # By removing every other frame
    def removeHalf(self,channels):
        train = []
        test = []
        for x in range(len(channels)):
            if x%2 :
                test.append(channels[x])
            else :
                train.append(channels[x])
        return train,test

    # Test the quality of the model
    def myTest(self,channels,oldLearner):
        trainingChannels,testChannels = self.removeHalf(channels)
        newLearner = lL(testChannels)
        distance = 0
        difference = 0
        for i in range(len(trainingChannels)-1):
            #some magical thing to work out where the point should be

            [x,y] = oldLearner.getPoint(i*2+1,0,1)

            guess = newLearner.getPose(array([x,y]),[0,1])
            xd = guess[0]-trainingChannels[i][0]
            yd = guess[1]-trainingChannels[i][1]
            zd = guess[2]-trainingChannels[i][2]
            distance+= sqrt( (xd*xd) + (yd*yd) + (zd*zd) )
            angleDiff = 0
            for j in range(3,len(trainingChannels[i])):
                angleDiff += abs(trainingChannels[i][j] - guess[j])

            angleDiff /= len(trainingChannels[i])-3
            difference += angleDiff

        distScore = (distance/(len(trainingChannels)-1))
        angleScore = (difference/(len(trainingChannels)-1))

        print distScore
        print angleScore

if __name__ == "__main__":
    r = runner()
    r.run('files/Swagger.bvh')














