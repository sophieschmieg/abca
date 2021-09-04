package fields.interfaces;

import java.util.Comparator;

public interface WellOrder<T extends Element<T>> extends MathSet<T>, Comparator<T> {
	public T least();
	public int compare(T t1, T t2);
	default public T max(T t1, T t2) {
		if (compare(t1, t2) < 0) { 
			return t2;
		}
		return t1;
	}
	default public T min(T t1, T t2) {
		if (compare(t1, t2) < 0) { 
			return t1;
		}
		return t2;
	}
}
