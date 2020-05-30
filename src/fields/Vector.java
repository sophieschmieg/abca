package fields;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.List;

public class Vector<T extends Element> implements Element {
	private T[] coeffs;

	@SuppressWarnings("unchecked")
	public Vector(List<T> coeffs) {
		this.coeffs = (T[]) Array.newInstance(coeffs.get(0).getClass(), coeffs.size());
		for (int i = 0; i < coeffs.size(); i++) {
			this.coeffs[i] = coeffs.get(i);
		}
	}

	@SafeVarargs
	public Vector(T... coeffs) {
		this.coeffs = Arrays.copyOf(coeffs, coeffs.length);
	}

	public int dimension() {
		return coeffs.length;
	}

	public T get(int index) {
		return coeffs[index - 1];
	}
	
	public List<T> asList() {
		return Arrays.asList(coeffs);
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Vector)) {
			return false;
		}
		@SuppressWarnings("unchecked")
		Vector<T> t = (Vector<T>) o;
		return Arrays.deepEquals(coeffs, t.coeffs);
	}

	@Override
	public int compareTo(Element o) {
		if (!(o instanceof Vector)) {
			throw new ArithmeticException();
		}
		@SuppressWarnings("unchecked")
		Vector<T> t = (Vector<T>) o;
		if (t.dimension() != this.dimension()) {
			return t.dimension() - this.dimension();
		}
		for (int i = 0; i < coeffs.length; i++) {
			int cmp = this.coeffs[i].compareTo(t.coeffs[i]);
			if (cmp != 0) {
				return cmp;
			}
		}
		return 0;
	}

	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append("(");
		buf.append(coeffs[0].toString());
		for (int i = 1; i < this.dimension(); i++) {
			buf.append(", ");
			buf.append(coeffs[i].toString());
		}
		buf.append(")");
		return buf.toString();
	}
}