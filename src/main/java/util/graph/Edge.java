package util.graph;

import fields.interfaces.Element;

public class Edge<V extends Element<V>, E extends Element<E>> implements Element<Edge<V, E>> {
	private boolean directed;
	private V firstVertex;
	private V secondVertex;
	private E weight;

	Edge(boolean directed, V firstVertex, V secondVertex, E weight) {
		this.directed = directed;
		this.firstVertex = firstVertex;
		this.secondVertex = secondVertex;
		this.weight = weight;
	}

	@Override
	public String toString() {
		return firstVertex + " --" + weight + "--" + (directed ? "> " : " ") + secondVertex;
	}

	public boolean directed() {
		return directed;
	}

	public V firstVertex() {
		return firstVertex;
	}

	public V secondVertex() {
		return secondVertex;
	}

	public V firstVertex(boolean direction) {
		if (direction) {
			return firstVertex;
		} else {
			return secondVertex;
		}
	}

	public V secondVertex(boolean direction) {
		if (direction) {
			return secondVertex;
		} else {
			return firstVertex;
		}
	}

	public boolean direction(V firstVertex) {
		if (this.firstVertex.equals(firstVertex)) {
			return true;
		} else if (this.secondVertex.equals(firstVertex)) {
			return false;
		}
		throw new ArithmeticException("Not a vertex");
	}

	public V otherVertex(V vertex) {
		if (direction(vertex)) {
			return secondVertex;
		}
		return firstVertex;
	}

	public E weight() {
		return weight;
	}

	@Override
	public int compareTo(Edge<V, E> o) {
		return hashCode() - o.hashCode();
	}
}
