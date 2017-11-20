/******************************************************************************
 * AUTHORS        : ____ & Jason Zhang
 * ASSIGNMENT #11 : BFS & DFS Directed Edges
 * CLASS          : CS1D
 * SECTION        : MW: 3:00pm
 * DUE DATE       : 11/20/2017
 *****************************************************************************/

#include "Graph.h"

/******************************************************************************
 * DFS & BFS
 * ----------------------------------------------------------------------------
 * This program utilizes a user defined Graph ADT in order to demonstrate a DFS
 * and a BFS on the vertices in a directed graph. The DFS will begin at Dallas
 * and visit all other vertices in the graph by choosing the shortest edge.
 * After the search is performed, the discovery, back, cross, and forward edges
 * will be displayed, as well as the total distance traveled on the discovery
 * edges.
 * Then, the BFS will begin at Dallas and visit all other vertices in the graph.
 * Afterwards, the discovery, back, and cross edges will be displayed, as well
 * as the total distance traveled.
 * ----------------------------------------------------------------------------
 * INPUT:
 *    <There is no input.>
 *
 * OUTPUT:
 *    DFS - Discovery Edges, Back Edges, Cross Edges, Forward Edges, and
 *          Total Distance Traveled
 *    BFS - Discovery Edges, Back Edges, Cross Edges, and Total Distance
 *          Traveled
 *
 *****************************************************************************/
int main()
{
	cout << "***************************************************************\n"
		 << "* NAME           : ____ & Jason Zhang\n"
		 << "* ASSIGNMENT #11 : DFS & BFS Directed Edges\n"
		 << "* CLASS          : CS1D - MW: 3:00pm\n"
		 << "* DUE DATE       : 11/20/2017\n*\n"
		 << "* DESCRIPTION: This program will demonstrate a DFS and BFS of a\n"
		 << "* directed graph using a user-defined Graph class. After each\n"
		 << "* traversal, the program will display the order that the vertices\n"
		 << "* were visited in, the total distance traveled, the discovery edges,\n"
		 << "* the back edges, the cross edges, and the forward edges (for DFS\n"
		 << "* only).\n"
		 <<	"***************************************************************\n\n";

	Graph graph; // Graph object.

	cout << "Adding vertices to the graph...\n";

	// Initializes the graph by reading from Map.txt.
	graph.initializeGraph();

	cout << "\n**********\n"
			"* PART A *\n"
			"**********\n";

	// Vector of city names that will hold the cities visited during the DFS, in
	// the order they were visited.
	vector<string> dfs;
	vector<string> ancestors;

	cout << "\nPerforming a DFS starting at Dallas:\n";

	// Performs a DFS on the graph starting at Dallas and stores the total
	// distance traveled.
	int dfsDistance = graph.DFS("Dallas", dfs, ancestors);

	for(unsigned int i = 0; i < dfs.size(); i++)
	{
		cout << dfs.at(i) << endl;
	}

	cout << "\nTotal Distance Traveled: " << dfsDistance << endl;

	// Vectors containing the discovery and back edges of the graph.
	vector<string> dfsDiscEdges;
    vector<string> dfsBackEdges;
    vector<string> dfsCrossEdges;
    vector<string> dfsForwardEdges;
    graph.dfsLabelEdges(dfs, dfsBackEdges,dfsDiscEdges,dfsCrossEdges,dfsForwardEdges);

    cout << "\nPrinting DFS discovery edges: \n";
    for(unsigned int i = 0; i < dfsDiscEdges.size(); i++)
        cout << dfsDiscEdges.at(i) << endl;

    cout << "\nPrinting DFS back edges: \n";
    for(unsigned int i = 0; i < dfsBackEdges.size(); i++)
        cout << dfsBackEdges.at(i) << endl;

    cout << "\nPrinting DFS cross edges: \n";
    for(unsigned int i = 0; i < dfsCrossEdges.size(); i++)
        cout << dfsCrossEdges.at(i) << endl;

    cout << "\nPrinting DFS forward edges: \n";
    for(unsigned int i = 0; i < dfsForwardEdges.size(); i++)
        cout << dfsForwardEdges.at(i) << endl;



	cout << "\n**********\n"
			"* PART B *\n"
			"**********\n\n";

    vector<string> bfs;

    cout << "Performing a BFS starting at Dallas: \n";
    int bfsDistance = graph.BFS("Dallas",bfs);

    for(unsigned int i = 0; i < bfs.size(); i++)
    {
        cout << bfs.at(i) << endl;
    }

    cout << "\nTotal Distance Traveled: " << bfsDistance << endl;

    vector<string> bfsDiscEdges;
    vector<string> bfsBackEdges;
    vector<string> bfsCrossEdges;
    graph.bfsLabelEdges(bfsBackEdges,bfsDiscEdges,bfsCrossEdges);

    cout << "\nPrinting BFS discovery edges: \n";
    for(unsigned int i = 0; i < bfsDiscEdges.size(); i++)
        cout << bfsDiscEdges.at(i) << endl;

    cout << "\nPrinting BFS back edges: \n";
    for(unsigned int i = 0; i < bfsBackEdges.size(); i++)
        cout << bfsBackEdges.at(i) << endl;

    cout << "\nPrinting BFS cross edges: \n";
    for(unsigned int i = 0; i < bfsCrossEdges.size(); i++)
        cout << bfsCrossEdges.at(i) << endl;

	return 0;
}
