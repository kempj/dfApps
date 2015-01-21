void init() {
    
}

int main(int argc; char* argv[]){
    //TODO: check to make sure numBlocks > numNodes
    int numBlocks, size;
    init();
    for(int i = 0; i < numBlocks; i++) {
        calcDiag(i,i);
        for(int j = i; k < numBlocksl j++) {
            calcCol(i,j);//split up over nodes (modulo or groups?)
        }
        //async send one block from each node that does the calculations to every node
        myNumRows = numBlocks/numNodes;
        for(int j = myNodeNum *myNumRows; j < (myNodeNum + 1) *myNumRows; j++) {
            calcRow(j);//split up over nodes (modulo)
        }
        int total = 0;
        while(total < myNumRows * (numBlocks - i)) {
            for(int j = i; j < myNumRows; j++) {
                for(int k = i; k < numBlocks; k++) {
                    if(received(j,k)) {
                        calcInner(j,k);
                        total++;
                    }
                }
            }
        }
    }
    return 0;
}
