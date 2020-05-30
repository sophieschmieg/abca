package util;

public class Pair<T, U> {
	private T first;
	private U second;
	
	public Pair(T first, U second) {
		super();
		this.first = first;
		this.second = second;
	}

	public T getFirst() {
		return first;
	}
	
	public U getSecond() {
		return second;
	}
	
	public String toString() {
		return "(" + first.toString() + ", " + second.toString() + ")";
	}
}
