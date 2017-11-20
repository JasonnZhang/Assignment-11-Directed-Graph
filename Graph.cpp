#include "Graph.h"

Graph::Graph()
{
	dfsDistance = 0;
}

Graph::~Graph()
{

}

bool Graph::empty()
{
	return graph.empty();
}

int Graph::size()
{
	return graph.size();
}

void Graph::initializeGraph()
{
	fstream fin; // File stream variable.

	// Opens the file.
	fin.open("Map.txt");

	// Reads in the contents of the file.
	while(!fin.eof())
	{
		string u; // Starting city.
		string v; // Ending city.
		int weight; // Distance between the cities.

		getline(fin, u, ',');
		getline(fin, v, ',');
		fin >> weight;
		fin.ignore(1000, '\n');

		// Inserts the edge to the graph.
		insertEdge(u, v, weight);
	}

	fin.close();
}

unsigned int Graph::findVertex(string city)
{
	unsigned int index = 0;
	bool found = false;

	// Searches for the city in the graph and returns the index of the graph it
	// was found, else returns the size of the graph.
	while(index < graph.size() && !found)
	{
		if(graph.at(index).city == city)
		{
			found = true;
		}
		else
		{
			index++;
		}
	}

	return index;
}

void Graph::insertVertex(string city)
{
	// Adds the vertex to the graph if it does not yet exist.
	if(findVertex(city) == graph.size())
	{
		Vertex newVertex;

		newVertex.city = city;
		newVertex.visited = false;

		cout << newVertex.city << endl;

		graph.push_back(newVertex);
	}
}

void Graph::insertEdge(string u, string v, int weight)
{
	// Returns the index of the vertex.
	unsigned int index = findVertex(u);

	// Adds the vertex to the graph if it does not yet exist and performs a
	// recursive call, else adds the edge to the edgeList of the specified vertex.
	if(index == graph.size())
	{
		insertVertex(u);

		insertEdge(u, v, weight);
	}
	else
	{
		Edge newEdge;

		newEdge.u = u;
		newEdge.v = v;
		newEdge.weight = weight;

		graph.at(index).edgeList.push_back(newEdge);
	}
}

vector<string> Graph::vertices()
{
	vector<string> cityNames; // Vector of city names.

	// Adds the names of the cities in the graph to the vector.
	for(unsigned int i = 0; i < graph.size(); i++)
	{
		cityNames.push_back(graph[i].city);
	}

	return cityNames;
}

vector<string> Graph::edges()
{
	vector<string> edgeList; // Vector of edges.

	// Adds the edges in the graph to the vector if the node has edges.
	for(unsigned int i = 0; i < graph.size(); i++)
	{
		if(graph.at(i).edgeList.size() != 0)
		{
			// Pushes the edge pair onto the vector in (u, v) format.
			for(unsigned int j = 0; j < graph.at(i).edgeList.size(); j++)
			{
				edgeList.push_back("(" + graph.at(i).edgeList.at(j).u + ", "
									   + graph.at(i).edgeList.at(j).v + ")");
			}
		}
	}

	return edgeList;
}

int Graph::DFS(string startingCity, vector<string> &dfs, vector<string> &ancestors)
{
	// Gets the graph index of the vertex being visited.
	int currVertex = findVertex(startingCity);

	// Visits the vertex.
	graph.at(currVertex).visited = true;

	// Adds the starting city to the vector
    dfs.push_back(startingCity);

    // Add starting city to ancestors vector
    ancestors.push_back(startingCity);

	// Performs a recursive call on itself to visit all vertices in the graph.
	if(verticesVisited() != graph.size())
	{
		// Gets the shortest outgoing edge
		Edge * nextEdge = smallestEdgeDFS(currVertex, dfs, ancestors);

		// Performs recursive call to visit the next closest city.
		DFS(nextEdge->v, dfs, ancestors);
	}

	return dfsDistance;
}

void Graph::dfsLabelEdges(vector<string> &dfs, vector<string> &backEdges, vector<string> &discEdges, vector<string> &crossEdges, vector<string> &forwardEdges)
{
    Edge * currEdge;
    string currString;

	// Adds the edges in the graph to the vector if the node has edges.
	for(unsigned int i = 0; i < graph.size(); i++)
	{
		int dfsIndex = findVertex(dfs.at(i));
		for(unsigned int j = 0; j < graph.at(dfsIndex).edgeList.size(); j++)
		{
			currEdge = &graph.at(dfsIndex).edgeList.at(j);
			currString = "(" + currEdge->u + ", " + currEdge->v + ")";
			switch (currEdge->dfsType)
            {
            case 0:
                backEdges.push_back(currString);
                break;
            case 1:
                discEdges.push_back(currString);
                break;
            case 2:
                // somehow figure out if it's a forward or cross edge
                forwardEdges.push_back(currString);
                break;
			}
		}
	}
}

void Graph::bfsLabelEdges(vector<string> &backEdges, vector<string> &discEdges, vector<string> &crossEdges)
{
    Edge * currEdge;
    string currString;

	// Adds the edges in the graph to the vector if the node has edges.
	for(unsigned int i = 0; i < graph.size(); i++)
	{
		for (unsigned int j=0; j<graph.at(i).edgeList.size(); j++)
		{
			currEdge = &graph.at(i).edgeList.at(j);
			currString = "(" + currEdge->u + ", " + currEdge->v + ")";
			if (currEdge->bfsType != 1)
            {
                if(graph.at(findVertex(currEdge->u)).bfsLevel > graph.at(findVertex(currEdge->v)).bfsLevel) {
                    currEdge->bfsType = 0;
                    backEdges.push_back(currString);
                }
                else {
                    currEdge->bfsType = 2;
                    crossEdges.push_back(currString);
                }
			}
			else {
                discEdges.push_back(currString);
			}
		}
	}
}

Edge * Graph::smallestEdgeDFS(int currVertex, vector<string> &dfs, vector<string> &ancestors)
{
	// Gets the size of the edgeList for the current vertex.
    int size = graph.at(currVertex).edgeList.size();

	// Searches for the next closest edge if all edges of the current vertex have
	// not been visited yet, else backtracks to find a vertex whose edges have
	// not all been discovered.

	if(edgesDiscovered(currVertex) != size)
	{
		// Edge list index of the closest city.
		int smallestIndex = 0;

		// Edge list index of the city whose distance is being compared to the
		// city at edgeList.at(smallestIndex).
		int compIndex = smallestIndex + 1;

		// If the other vertex in the edgeList.at(0) is not visited and not an
		// ancestor, mark the edge as a cross edge (temporarily).
		// Note: Index 0 needs to be checked separately here, since compIndex
		// will start at index 1
        if (!graph.at(findVertex(graph.at(currVertex).edgeList.at(0).v)).visited) {
            vector<string>::iterator ancestorsIt =
                find(ancestors.begin(),ancestors.end(),graph.at(currVertex).edgeList.at(0).v);
            if (ancestorsIt == ancestors.end())
                graph.at(currVertex).edgeList.at(0).dfsType = 2;
        }

		// Finds the next closest city that has not been visited yet.
		while(compIndex < size)
		{
			// Gets the graph index of the next closest city.
			int smallestVertex = findVertex(graph.at(currVertex).edgeList.at(smallestIndex).v);

			// Gets the graph index of the city in the edge list being comapred
			// to the city at edgeList.at(smallestIndex).
			int compVertex = findVertex(graph.at(currVertex).edgeList.at(compIndex).v);


			// If the vertex at graph.at(smallestVertex) has already been visited,
			// increments smallest index and does nothing, else checks if the
			// vertex it is being compared to has been visited.
			if(graph.at(smallestVertex).visited)
			{
				smallestIndex++;
			}
			else
			{
				// If the vertex smallestVertex is being compared to has not been
				// visited, compares their weights, else does nothing.
				if(!(graph.at(compVertex).visited))
				{
					if(graph.at(currVertex).edgeList.at(smallestIndex).weight >=
					   graph.at(currVertex).edgeList.at(compIndex).weight)
					{
						smallestIndex = compIndex;
					}
				}

				// If the other vertex in the edge is not visited and not an ancestor,
                // mark the edge as a cross edge (temporarily).
                vector<string>::iterator ancestorsIt =
                    find(ancestors.begin(),ancestors.end(),graph.at(currVertex).edgeList.at(compIndex).v);
                if (ancestorsIt == ancestors.end())
                    graph.at(currVertex).edgeList.at(compIndex).dfsType = 2;
			}

			// Increments compIndex so that it is always the after smallestIndex.
			compIndex++;
		}

		// Marks the edge that has the smallest weight as a discovery edge.
		graph.at(currVertex).edgeList.at(smallestIndex).dfsType = 1;

		// Adds the distance to the overall distance traveled.
		dfsDistance += graph.at(currVertex).edgeList.at(smallestIndex).weight;

		return &graph.at(currVertex).edgeList.at(smallestIndex);
	}
	else
	{
		// Iterator that gets the location of the current city in the vector of
		// names that contains the cities in the order they were visited.
		vector<string>::iterator dfsIt = find(dfs.begin(), dfs.end(),
										 	  graph.at(currVertex).city);

		// Decrements the iterator to the previous city visited.
		dfsIt--;

		// Finds the graph index of the previous city visited.
		int backIndex = findVertex(*dfsIt);

		// Preforms a recursive call to check if the previous city visited has
		// any unvisited edges to continue the DFS.
		ancestors.pop_back();
		return smallestEdgeDFS(backIndex, dfs, ancestors);
	}
}

int Graph::BFS(string startingCity, vector<string> &bfs)
{
    // Reset the graph, this should be its own function
    for (unsigned int i=0; i<graph.size(); i++) {
        graph.at(i).visited = false;
        graph.at(i).bfsLevel = 0;
        for (unsigned int j=0; j< graph.at(i).edgeList.size(); j++)
            graph.at(i).edgeList.at(j).bfsType = 1;
    }

    // Get the graph index of the vertex being visited.
	int currVertex = findVertex(startingCity);

	// Visit the starting vertex
	graph.at(currVertex).visited = true;
	graph.at(currVertex).bfsLevel = 0;

	// Create a newLevel vector to hold the current level, containing only
	// the 1st vertex, and add the first vertex to the bfs vector
	vector<int> newLevel;
	newLevel.push_back(currVertex);
	bfs.push_back(startingCity);

	// Start recursion
	return BFSRecur(bfs, newLevel, 0);
}

int Graph::BFSRecur(vector<string> &bfs, vector<int> previousLevel, int bfsLevel)
{
    vector<int> newLevel;
    vector<int> currLevel;
    int levelDistance = 0;

    vector<Edge> * currEdgeList;
    Vertex * startingVertex;
    Vertex * currVertex;
    int currVertexID;

    // Iterate through the previous level
    for (unsigned int i=0; i<previousLevel.size(); i++) {
        startingVertex = &graph.at(previousLevel.at(i));
        currEdgeList = &startingVertex->edgeList;

        // Iterate through the ith vertex's edge list
        for (unsigned int j=0; j<currEdgeList->size(); j++)
        {
            // Add all non-visited levels to the next level, in closest order
            currVertexID = findVertex(currEdgeList->at(j).v);
            currVertex = &graph.at(currVertexID);
            if (!currVertex->visited)
            {
                // Add the edge length to the distance, including return trip
                levelDistance += currEdgeList->at(j).weight;

                // Mark the vertex as visited and the edge as a discovery edge
                currVertex->visited = true;
                currVertex->bfsLevel = bfsLevel;
                currEdgeList->at(j).bfsType = 1;

                // Insert the current vertex in the sorted position
                bool inserted = false; // could do the same thing by changing the visited variable, but this is clearer
                for (unsigned int k=0; k<currLevel.size() && !inserted; k++)
                {
                    if (currEdgeList->at(j).weight < distance(startingVertex, &graph.at(currLevel.at(k))))
                    {
                        currLevel.insert(currLevel.begin()+k,currVertexID);
                        inserted = true;
                    }
                }
                if (!inserted)
                    currLevel.push_back(findVertex(currVertex->city));
            }
            else {
                currEdgeList->at(j).bfsType = 0;
            }
        }
        // Add the current level vertices to the the end of the bfs vector
        for (unsigned int m=0; m<currLevel.size(); m++) {
            bfs.push_back(graph.at(currLevel.at(m)).city);
        }

        // Add the current level vertices to the end of the newLevel vector
        newLevel.insert(newLevel.end(),currLevel.begin(),currLevel.end());
        currLevel.clear();
    }

    // If still has vertices, do recursive call
    if (newLevel.size() > 0)
        return levelDistance + BFSRecur(bfs, newLevel, bfsLevel+1);
    else
        return levelDistance;
}

int Graph::distance(Vertex * v1, Vertex * v2)
{
    // find connecting edge
    for (unsigned int i=0; i<v1->edgeList.size(); i++) {
        if (v1->edgeList.at(i).v == v2->city)
            return v1->edgeList.at(i).weight;
    }
    return -1;
}

string Graph::otherVertex(Edge currEdge, string startingCity)
{
    if(currEdge.u == startingCity)
        return currEdge.v;
    else
        return currEdge.u;
}

unsigned int Graph::verticesVisited()
{
	int numVisited = 0; // Number of vertices visited.

	// Iterator to the first vertex in the graph.
	vector<Vertex>::iterator graphIt = graph.begin();

	// Counts the number of vertices in the graph that are marked as visited.
	while(graphIt != graph.end())
	{
		if(graphIt->visited)
		{
			numVisited++;
		}

		graphIt++;
	}

	return numVisited;
}

unsigned int Graph::edgesDiscovered(int currVertex)
{
	int numVisited = 0; // Number of edges discovered.

	// Counts the number of edges at the current vertex that have been discovered.
	for(unsigned int i = 0; i < graph.at(currVertex).edgeList.size(); i++)
	{
		if(graph.at(findVertex(graph.at(currVertex).edgeList.at(i).v)).visited)
		{
			numVisited++;
		}
	}

	return numVisited;
}
