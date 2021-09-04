package fields.interfaces;

public interface Element<T extends Element<T>> extends Comparable<T> {
	public boolean equals(Object o);
}
