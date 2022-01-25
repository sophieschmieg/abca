package util.graph;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.interfaces.Element;
import fields.interfaces.Group;
import fields.interfaces.MathMap;
import fields.interfaces.Monoid;
import fields.vectors.Matrix;
import util.ConstantMap;
import util.Identity;

public class Graph<V extends Element<V>, E extends Element<E>> {
	private Set<Edge<V, E>> edgeSet;
	private SortedMap<V, Set<Edge<V, E>>> adjacencyLists;
	private boolean directed;
	private boolean pathConnected;
	private boolean pathConnectionComputed;
	private boolean multiple;
	private boolean multipleComputed;
	private List<Graph<V, E>> connectedComponents;
	private Matrix<E> adjacencyMatrix;
	private Map<V, Map<V, Optional<Path>>> shortestPaths;
	private Map<V, Map<V, Optional<Path>>> paths;
	private List<Cycle> simpleCycles;
	private Graph<V, E> asUndirectedGraph;

	private Graph(boolean directed, Set<Edge<V, E>> edgeSet, SortedMap<V, Set<Edge<V, E>>> adjacencyLists) {
		this.directed = directed;
		this.edgeSet = edgeSet;
		this.adjacencyLists = adjacencyLists;
		this.multipleComputed = false;
	}

	public boolean isDirected() {
		return directed;
	}

	public Set<V> vertexSet() {
		return adjacencyLists.keySet();
	}

	public Set<Edge<V, E>> edgeSet() {
		return edgeSet;
	}

	public Set<Edge<V, E>> edges(V vertex) {
		return adjacencyLists.get(vertex);
	}

	public Graph<V, E> asUndirectedGraph() {
		if (!directed) {
			return this;
		}
		if (asUndirectedGraph == null) {
			GraphBuilder<V, E> builder = new GraphBuilder<>(false);
			builder.addAll(this);
			asUndirectedGraph = builder.build();
		}
		return asUndirectedGraph;
	}

	private class UnionFindEntry {
		private UnionFindEntry parent;
		private int weight;
		private V vertex;

		private UnionFindEntry(V vertex) {
			this.vertex = vertex;
			this.weight = 1;
			this.parent = null;
		}

		private UnionFindEntry find() {
			if (parent == null) {
				return this;
			}
			parent = parent.find();
			return parent;
		}

		private void merge(UnionFindEntry with) {
			UnionFindEntry thisParent = find();
			UnionFindEntry otherParent = with.find();
			if (thisParent == otherParent) {
				return;
			}
			if (thisParent.weight < otherParent.weight) {
				thisParent.parent = otherParent;
				otherParent.weight += thisParent.weight;
			} else {
				otherParent.parent = thisParent;
				thisParent.weight += otherParent.weight;
			}
		}
	}

	public List<Graph<V, E>> connectedComponents() {
		if (connectedComponents == null) {
			Map<V, UnionFindEntry> unionFind = new TreeMap<>();
			for (V vertex : vertexSet()) {
				unionFind.put(vertex, new UnionFindEntry(vertex));
			}
			for (Edge<V, E> edge : edgeSet()) {
				unionFind.get(edge.firstVertex()).merge(unionFind.get(edge.secondVertex()));
			}
			Map<V, GraphBuilder<V, E>> builders = new TreeMap<>();
			for (UnionFindEntry entry : unionFind.values()) {
				builders.putIfAbsent(entry.find().vertex, new GraphBuilder<>(directed));
				builders.get(entry.find().vertex).addVertex(entry.vertex);
			}
			for (Edge<V, E> edge : edgeSet()) {
				builders.get(unionFind.get(edge.firstVertex()).find().vertex).addEdge(edge.firstVertex(),
						edge.secondVertex(), edge.weight());
			}
			connectedComponents = new ArrayList<>();
			for (GraphBuilder<V, E> builder : builders.values()) {
				Graph<V, E> component = builder.build();
				component.connectedComponents = Collections.singletonList(component);
				connectedComponents.add(component);
			}
		}
		return connectedComponents;
	}

	public boolean isConnected() {
		return connectedComponents().size() == 1;
	}

	public boolean isPathConnected() {
		if (!pathConnectionComputed) {
			if (!directed) {
				pathConnected = isConnected();
			} else {
				if (!isConnected()) {
					pathConnected = false;
				} else {
					pathConnected = true;
					fromLoop: for (V from : vertexSet()) {
						for (V to : vertexSet()) {
							Optional<Path> path = findPath(from, to);
							if (path.isEmpty()) {
								pathConnected = false;
								break fromLoop;
							}
						}
					}
				}
			}
			pathConnectionComputed = true;
		}
		return pathConnected;
	}

	public int cyclomaticNumber() {
		return edgeSet().size() - vertexSet().size() + connectedComponents().size();
	}
	
	private <T extends Element<T>>Path reconstructPath(Dijkstra<T> end) {
		Dijkstra<T> current = end;
		List<Edge<V, E>> edges = new ArrayList<>();
		List<Boolean> directions = new ArrayList<>();
		while (current.prev != null) {
			edges.add(current.edge);
			directions.add(current.direction);
			current = current.prev;
		}
		Collections.reverse(edges);
		Collections.reverse(directions);
		return new Path(edges, directions);
	}

	private void computeSimpleCyclesDirected() {
		if (simpleCycles != null) {
			return;
		}
		simpleCycles = new ArrayList<>();
		Set<Cycle> knownCycles = new TreeSet<>();
		for (V vertex : vertexSet()) {
			Deque<Dijkstra<E>> stack = new LinkedList<>();
			Set<V> visited = new TreeSet<>();
			stack.push(new Dijkstra<>(vertex, null));
			while (!stack.isEmpty()) {
				Dijkstra<E> current = stack.pop();
				if (current.prev != null && current.vertex.equals(vertex)) {
					Cycle cycle = reconstructPath(current).asCycle();
					if (!cycle.isSimple()) {
						continue;
					}
					if (knownCycles.contains(cycle)) {
						continue;
					}
					simpleCycles.add(cycle);
					knownCycles.add(cycle);
				}
				if (visited.contains(current.vertex)) {
					continue;
				}
				visited.add(current.vertex);
				for (Edge<V, E> edge : edges(current.vertex)) {
					stack.push(new Dijkstra<>(edge.secondVertex(), edge, true, null, current));
				}
			}
		}
	}

	private void computeSimpleCyclesUndirected() {
		if (simpleCycles != null) {
			return;
		}
		simpleCycles = new ArrayList<>();
		for (Graph<V, E> component : connectedComponents()) {
			Set<Edge<V, E>> seen = new HashSet<>();
			int cycles = component.cyclomaticNumber();
			V start = component.vertexSet().iterator().next();
			Map<V, Dijkstra<E>> paths = new TreeMap<>();
			Deque<Dijkstra<E>> queue = new LinkedList<>();
			queue.add(new Dijkstra<E>(start, null));
			int cyclesFound = 0;
			while (!queue.isEmpty() && cyclesFound < cycles) {
				Dijkstra<E> current = queue.poll();
				if (paths.containsKey(current.vertex)) {
					Dijkstra<E> old = paths.get(current.vertex);
					cyclesFound++;
					Path currentToStart = reconstructPath(current);
					Path oldToStart = reconstructPath(old);
					Cycle cycle = simplifyCycle(concat(currentToStart, reversePath(oldToStart)).asCycle());
					simpleCycles.add(cycle);
					continue;
				}
				paths.put(current.vertex, current);
				for (Edge<V, E> edge : component.edges(current.vertex)) {
					if (seen.contains(edge)) {
						continue;
					}
					seen.add(edge);
					boolean direction = edge.direction(current.vertex);
					queue.add(new Dijkstra<>(edge.secondVertex(direction), edge, direction, null, current));
				}
			}
			if (cyclesFound != cycles) {
				throw new ArithmeticException("Could not find enough cycles!");
			}
		}
	}

	public List<Cycle> simpleCycles() {
		if (simpleCycles == null) {
			if (directed) {
				computeSimpleCyclesDirected();
			} else {
				computeSimpleCyclesUndirected();
			}
		}
		return simpleCycles;
	}

	public Group<Cycle> cycleGroup() {
		return new Group<>() {
			private Cycle emptyCycle = new Cycle(new Path(Collections.emptyList(), Collections.emptyList()));

			@Override
			public Exactness exactness() {
				return Exactness.EXACT;
			}

			@Override
			public Cycle getRandomElement() {
				List<Cycle> simple = simpleCycles();
				return simple.get(new Random().nextInt(simple.size()));
			}

			@Override
			public boolean isFinite() {
				return false;
			}

			@Override
			public BigInteger getNumberOfElements() throws InfinityException {
				throw new InfinityException();
			}

			@Override
			public Iterator<Cycle> iterator() {
				throw new InfinityException();
			}

			@Override
			public Cycle neutral() {
				return emptyCycle;
			}

			@Override
			public Cycle inverse(Cycle t) {
				return reversePath(t.asPath()).asCycle();
			}

			@Override
			public Cycle operate(Cycle t1, Cycle t2) {
				if (t1.empty()) {
					return t2;
				}
				if (t2.empty()) {
					return t1;
				}
				Optional<Path> connection = findPath(t1.asPath().last(), t2.asPath().first());
				if (connection.isEmpty()) {
					throw new ArithmeticException("Not a connected graph!");
				}
				Path path = concat(t1.asPath(), connection.get());
				path = concat(path, t2.asPath());
				path = concat(path, reversePath(connection.get()));
				return simplifyCycle(path.asCycle());
			}
		};
	}

	private class Dijkstra<T extends Element<T>> implements Comparable<Dijkstra<T>> {
		private V vertex;
		private T totalWeight;
		private Edge<V, E> edge;
		private boolean direction;
		private Dijkstra<T> prev;

		private Dijkstra(V start, T zero) {
			this.vertex = start;
			this.totalWeight = zero;
		}

		private Dijkstra(V vertex, Edge<V, E> edge, boolean direction, T totalWeight, Dijkstra<T> prev) {
			this.vertex = vertex;
			this.edge = edge;
			this.direction = direction;
			this.totalWeight = totalWeight;
			this.prev = prev;
		}

		@Override
		public String toString() {
			return vertex + " prev " + edge;
		}

		@Override
		public int compareTo(Dijkstra<T> o) {
			if (totalWeight == null && o.totalWeight == null) {
				return 0;
			}
			if (totalWeight == null) {
				return 1;
			}
			if (o.totalWeight == null) {
				return -1;
			}
			return totalWeight.compareTo(o.totalWeight);
		}
	}

	public Optional<Path> findPath(V from, V to) {
		if (paths == null) {
			paths = findAllShortestPaths(Integers.z().getAdditiveGroup(), new ConstantMap<>(Integers.z().one()));
		}
		return paths.get(from).get(to);
	}

	public Optional<Path> findShortestPath(V from, V to, Monoid<E> orderedWeightMonoid) {
		if (shortestPaths == null) {
			shortestPaths = findAllShortestPaths(orderedWeightMonoid, new Identity<>());
		}
		return shortestPaths.get(from).get(to);
	}

	private <T extends Element<T>> Map<V, Dijkstra<T>> findAllShortestPathDijkstra(V from,
			Monoid<T> orderedWeightMonoid, MathMap<E, T> weightMap) {
		PriorityQueue<Dijkstra<T>> queue = new PriorityQueue<>(vertexSet().size());
		Map<V, Dijkstra<T>> shortestKnownPath = new TreeMap<>();
		queue.add(new Dijkstra<>(from, orderedWeightMonoid.neutral()));
		while (!queue.isEmpty()) {
			Dijkstra<T> current = queue.poll();
			if (shortestKnownPath.containsKey(current.vertex)
					&& shortestKnownPath.get(current.vertex).compareTo(current) <= 0) {
				continue;
			} else {
				shortestKnownPath.put(current.vertex, current);
			}
			for (Edge<V, E> edge : edges(current.vertex)) {
				V otherVertex = edge.firstVertex();
				boolean direction = false;
				if (current.vertex.equals(otherVertex)) {
					otherVertex = edge.secondVertex();
					direction = true;
				}
				queue.add(new Dijkstra<>(otherVertex, edge, direction,
						orderedWeightMonoid.operate(current.totalWeight, weightMap.evaluate(edge.weight())), current));
			}
		}
		return shortestKnownPath;
	}

	public <T extends Element<T>> Map<V, Map<V, Optional<Path>>> findAllShortestPaths(Monoid<T> orderedWeightMonoid,
			MathMap<E, T> weightMap) {
		Map<V, Map<V, Optional<Path>>> result = new TreeMap<>();
		for (V from : vertexSet()) {
			result.put(from, new TreeMap<>());
			Map<V, Dijkstra<T>> shortestPathsDijkstra = findAllShortestPathDijkstra(from, orderedWeightMonoid,
					weightMap);
			for (V to : vertexSet()) {
				if (!shortestPathsDijkstra.containsKey(to)) {
					result.get(from).put(to, Optional.empty());
				}
				Dijkstra<T> current = shortestPathsDijkstra.get(to);
				result.get(from).put(to, Optional.of(reconstructPath(current)));
			}
		}
		return result;
	}

	public <T extends Element<T>> Optional<Path> findShortestPath(V from, V to, Monoid<T> orderedWeightMonoid,
			MathMap<E, T> weightMap) {
		return findAllShortestPaths(orderedWeightMonoid, weightMap).get(from).get(to);
	}

	public Path concat(Path first, Path second) {
		if (first.empty()) {
			return second;
		}
		if (second.empty()) {
			return first;
		}
		if (!first.last().equals(second.first())) {
			throw new ArithmeticException("Paths are not concatable!");
		}
		List<Edge<V, E>> edges = new ArrayList<>();
		List<Boolean> directions = new ArrayList<>();
		edges.addAll(first.edges());
		edges.addAll(second.edges());
		directions.addAll(first.directions());
		directions.addAll(second.directions());
		return new Path(edges, directions);
	}

	public Cycle simplifyCycle(Cycle cycle) {
		Path simplifiedPath = simplifyPath(cycle.asPath());
		List<Edge<V, E>> edges = new ArrayList<>();
		List<Boolean> directions = new ArrayList<>();
		edges.addAll(simplifiedPath.edges());
		directions.addAll(simplifiedPath.directions());
		while (edges.get(0) == edges.get(edges.size() - 1)
				&& directions.get(0) != directions.get(directions.size() - 1)) {
			edges = edges.subList(1, edges.size() - 1);
			directions = directions.subList(1, directions.size() - 1);
		}
		return new Path(edges, directions).asCycle();
	}

	public Path simplifyPath(Path path) {
		List<Edge<V, E>> edges = new ArrayList<>();
		List<Boolean> directions = new ArrayList<>();
		for (int i = 0; i < path.edges().size(); i++) {
			if (edges.isEmpty()) {
				edges.add(path.edges().get(i));
				directions.add(path.directions().get(i));
				continue;
			}
			Edge<V, E> lastEdge = edges.get(edges.size() - 1);
			boolean lastDirection = directions.get(directions.size() - 1);
			if (lastEdge == path.edges().get(i) && lastDirection != path.directions().get(i)) {
				edges.remove(edges.size() - 1);
				directions.remove(directions.size() - 1);
			} else {
				edges.add(path.edges().get(i));
				directions.add(path.directions().get(i));
			}
		}
		return new Path(edges, directions);
	}

	public Path reversePath(Path path) {
		List<Edge<V, E>> edges = new ArrayList<>();
		List<Boolean> directions = new ArrayList<>();
		edges.addAll(path.edges());
		Collections.reverse(edges);
		for (int i = path.directions().size() - 1; i >= 0; i--) {
			directions.add(!path.directions().get(i));
		}
		return new Path(edges, directions);
	}

	public boolean isMultiple() {
		if (!multipleComputed) {
			for (V vertex : vertexSet()) {
				Set<V> vertexSeen = new TreeSet<>();
				for (Edge<V, E> edge : edges(vertex)) {
					V secondVertex;
					if (edge.firstVertex().equals(vertex)) {
						secondVertex = edge.secondVertex();
					} else {
						if (directed) {
							continue;
						}
						secondVertex = edge.firstVertex();
					}
					if (vertexSeen.contains(secondVertex)) {
						multiple = true;
						multipleComputed = true;
						return true;
					}
					vertexSeen.add(secondVertex);
				}
			}
			multiple = false;
			multipleComputed = true;
		}
		return multiple;
	}

	public Matrix<E> adjacencyMatrix(E zero) {
		if (adjacencyMatrix == null) {
			List<List<E>> edgeWeights = new ArrayList<>();
			int numVertex = vertexSet().size();
			Map<V, Integer> vertexIndeces = new TreeMap<>();
			int i = 0;
			for (V vertex : vertexSet()) {
				vertexIndeces.put(vertex, i);
				i++;
				List<E> row = new ArrayList<>();
				for (int j = 0; j < numVertex; j++) {
					row.add(zero);
				}
				edgeWeights.add(row);
			}
			for (Edge<V, E> edge : edgeSet) {
				i = vertexIndeces.get(edge.firstVertex());
				int j = vertexIndeces.get(edge.secondVertex());
				E weight = edgeWeights.get(i).get(j);
				if (!weight.equals(zero)) {
					throw new ArithmeticException("multiple graph!");
				}
				edgeWeights.get(i).set(j, edge.weight());
				if (!directed) {
					edgeWeights.get(j).set(i, edge.weight());
				}
			}
			adjacencyMatrix = new Matrix<>(edgeWeights);
		}
		return adjacencyMatrix;
	}

	public class Path extends AbstractElement<Path> {
		private List<V> path;
		private List<E> edgePath;
		private List<Edge<V, E>> edges;
		private List<Boolean> directions;

		private Path(List<Edge<V, E>> edges, List<Boolean> directions) {
			this.path = new ArrayList<>();
			this.edgePath = new ArrayList<>();
			this.edges = edges;
			this.directions = directions;
			if (edges.size() != directions.size()) {
				throw new ArithmeticException("path size mismatch!");
			}
			for (int i = 0; i < edges.size(); i++) {
				if (directed && !directions.get(i)) {
					throw new ArithmeticException("Path not following the direction!");
				}
				V firstVertex = edges.get(i).firstVertex(directions.get(i));
				V secondVertex = edges.get(i).secondVertex(directions.get(i));
				if (i == 0) {
					path.add(firstVertex);
				}
				path.add(secondVertex);
				if (!path.get(i).equals(firstVertex)) {
					throw new ArithmeticException("Path not connected!");
				}
				edgePath.add(edges.get(i).weight());
			}
		}

		@Override
		public String toString() {
			if (empty()) {
				return "e";
			}
			StringBuilder build = new StringBuilder();
			build.append(path.get(0));
			for (int i = 0; i < edgePath.size(); i++) {
				build.append(" -- " + edgePath.get(i) + " -- " + path.get(i + 1));
			}
			return build.toString();
		}

		public boolean empty() {
			return path.isEmpty();
		}

		public boolean isCycle() {
			return empty() || first().equals(last());
		}

		public Cycle asCycle() {
			return new Cycle(this);
		}

		public V first() {
			return path.get(0);
		}

		public V last() {
			return path.get(path.size() - 1);
		}

		public List<V> vertexPath() {
			return path;
		}

		public List<E> edgePath() {
			return edgePath;
		}

		public List<Edge<V, E>> edges() {
			return edges;
		}

		public List<Boolean> directions() {
			return directions;
		}

		@Override
		public int compareTo(Path o) {
			if (path.size() != o.path.size()) {
				return path.size() - o.path.size();
			}
			for (int i = 0; i < path.size(); i++) {
				int cmp = path.get(i).compareTo(o.path.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			for (int i = 0; i < path.size(); i++) {
				int cmp = edgePath.get(i).compareTo(o.edgePath.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			for (int i = 0; i < edgePath.size(); i++) {
				int cmp = edgePath.get(i).compareTo(o.edgePath.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			for (int i = 0; i < edges.size(); i++) {
				int cmp = edges.get(i).compareTo(o.edges.get(i));
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
	}

	public class Cycle extends AbstractElement<Cycle> {
		private Path path;

		private Cycle(Path path) {
			this.path = path;
			if (!path.empty() && !path.first().equals(path.last())) {
				throw new ArithmeticException("Not a cycle!");
			}
		}

		@Override
		public String toString() {
			return path.toString();
		}

		public boolean empty() {
			return path.empty();
		}

		public boolean isSimple() {
			if (empty()) {
				return false;
			}
			Set<V> seen = new TreeSet<>();
			for (V vertex : vertexPath()) {
				if (seen.contains(vertex)) {
					return false;
				}
				seen.add(vertex);
			}
			return true;
		}

		public Path asPath() {
			return path;
		}

		public List<V> vertexPath() {
			return path.vertexPath().subList(0, path.edges().size());
		}

		public List<E> edgePath() {
			return path.edgePath();
		}

		public List<Edge<V, E>> edges() {
			return path.edges();
		}

		public List<Boolean> directions() {
			return path.directions();
		}

		@Override
		public int compareTo(Cycle o) {
			if (path.path.size() != o.path.path.size()) {
				return path.path.size() - o.path.path.size();
			}
			if (empty()) {
				return 0;
			}
			Edge<V, E> firstEdge = path.edges.get(0);
			List<Integer> starts = new ArrayList<>();
			for (int i = 0; i < path.edges.size(); i++) {
				if (o.path.edges.get(i) == firstEdge) {
					starts.add(i);
				}
			}
			if (starts.isEmpty()) {
				return path.compareTo(o.path);
			}
			startLoop: for (int start : starts) {
				for (int i = 0; i < path.edges.size(); i++) {
					if (path.edges.get(i) != path.edges.get((i + start) % path.edges.size())) {
						continue startLoop;
					}
				}
				return 0;
			}
			return path.compareTo(o.path);
		}
	}

	public static class GraphBuilder<V extends Element<V>, E extends Element<E>> {
		private Set<Edge<V, E>> edgeSet;
		private SortedMap<V, Set<Edge<V, E>>> adjacencyLists;
		private boolean directed;
		private boolean built;

		public GraphBuilder(boolean directed) {
			this.directed = directed;
			this.edgeSet = new TreeSet<>();
			this.adjacencyLists = new TreeMap<>();
			built = false;
		}

		public void addVertex(V vertex) {
			adjacencyLists.putIfAbsent(vertex, new TreeSet<>());
		}

		public void addEdge(V firstVertex, V secondVertex, E weight) {
			addVertex(firstVertex);
			addVertex(secondVertex);
			Edge<V, E> newEdge = new Edge<>(directed, firstVertex, secondVertex, weight);
			edgeSet.add(newEdge);
			adjacencyLists.get(firstVertex).add(newEdge);
			if (!directed) {
				adjacencyLists.get(secondVertex).add(newEdge);
			}
		}

		public void addAll(Graph<V, E> graph) {
			for (V vertex : graph.vertexSet()) {
				addVertex(vertex);
			}
			for (Edge<V, E> edge : graph.edgeSet()) {
				addEdge(edge.firstVertex(), edge.secondVertex(), edge.weight());
			}
		}

		public Graph<V, E> build() {
			if (built) {
				throw new ArithmeticException("Graph already built!");
			}
			built = true;
			return new Graph<>(directed, edgeSet, adjacencyLists);
		}
	}
}
