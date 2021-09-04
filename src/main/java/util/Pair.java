package util;

public class Pair<T extends Comparable<? super T>, U extends Comparable<? super U>> implements Comparable<Pair<T, U>> {
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

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Pair<?, ?>)) {
			return false;
		}
		Pair<?, ?> other = (Pair<?, ?>) o;
		return this.first.equals(other.first) && this.second.equals(other.second);
	}

	@Override
	public int compareTo(Pair<T, U> o) {
		int cmp = first.compareTo(o.first);
		if (cmp != 0) {
			return cmp;
		}
		return second.compareTo(o.second);
	}
}
