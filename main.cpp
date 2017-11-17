/******************************************************************************
 * AUTHORS        : Jesus Sanchez & Jason Zhang
 * ASSIGNMENT #10 : BFS & DFS
 * CLASS          : CS1D
 * SECTION        : MW: 3:00pm
 * DUE DATE       : 11/13/2017
 *****************************************************************************/

#include "Graph.h"

/******************************************************************************
 * DFS & BFS
 * ----------------------------------------------------------------------------
 * This program utilizes a user defined Graph ADT in order to demonstrate a DFS
 * and a BFS on the vertices in the graph. The DFS will begin at Dallas and
 * visit all other vertices in the graph by choosing the edge with the smallest
 * mileage. After the search is performed, the discovery and back edges will be
 * displayed, as well as the total distance traveled on the discovery edges.
 * Then, the BFS will be performed and the discovery and cross edges will be
 * printed, as well as the total distance traveled.
 * ----------------------------------------------------------------------------
 * INPUT:
 *    <There is no input.>
 *
 * OUTPUT:
 *    DFS - Discovery Edges, Back Edges, and Total Distance Traveled
 *    BFS - Discovery Edges, Cross Edges, and Total Distance Traveled
 *
 *****************************************************************************/
int main()
{
	cout << "***************************************************************\n"
		 << "* NAME           : Jesus Sanchez & Jason Zhang\n"
		 << "* ASSIGNMENT #10 : DFS & BFS\n"
		 << "* CLASS          : CS1D - MW: 3:00pm\n"
		 << "* DUE DATE       : 11/13/2017\n*\n"
		 << "* DESCRIPTION: This program will demonstrate a DFS and BFS using\n"
		 << "* a user-defined Graph class. After each traversal, the order\n"
		 <<	"* that the vertices were visited in will be displayed, as well\n"
		 << "* as the discovery edges and back/cross edges.\n"
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
	vector<string> dfsDiscoveryEdges = graph.getDiscoveryEdges(dfs);
	vector<string> dfsBackEdges = graph.getBackEdges(dfs);

	cout << "\nPrinting DFS discovery edges:\n";
	for(unsigned int i = 0; i < dfsDiscoveryEdges.size(); i++)
	{
		cout << dfsDiscoveryEdges.at(i) << endl;
	}

	cout << "\nPrinting DFS back edges:\n";
	for(unsigned int i = 0; i < dfsBackEdges.size(); i++)
	{
		cout << dfsBackEdges.at(i) << endl;
	}

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

    vector<string> bfsBackEdges;
    vector<string> bfsCrossEdges;
    graph.bfsLabelEdges(bfsBackEdges,bfsCrossEdges);

    cout << "\nPrinting BFS back edges: \n";

    for(unsigned int i = 0; i < bfsBackEdges.size(); i++)
    {
        cout << bfsBackEdges.at(i) << endl;
    }

    cout << "\nPrinting BFS cross edges: \n";

    for(unsigned int i = 0; i < bfsCrossEdges.size(); i++)
    {
        cout << bfsCrossEdges.at(i) << endl;
    }

	return 0;
}
