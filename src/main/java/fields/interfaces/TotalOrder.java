package fields.interfaces;

import java.util.Comparator;

public interface TotalOrder<T extends Element<T>> extends MathSet<T>, Comparator<T> {

	@Override
	default int compare(T o1, T o2) {
		return o1.compareTo(o2);
	}

	default T max(T t1, T t2) {
		if (t1.compareTo(t2) < 0) {
			return t2;
		}
		return t1;
	}

	default T min(T t1, T t2) {
		if (t1.compareTo(t2) < 0) {
			return t1;
		}
		return t2;
	}

}
