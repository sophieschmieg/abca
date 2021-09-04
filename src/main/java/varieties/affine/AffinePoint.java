package varieties.affine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;

public class AffinePoint<T extends Element<T>> extends AbstractElement<AffinePoint<T>> {
	private int dim;
	private List<T> coords;

	public AffinePoint(Field<T> field, T coord) {
		this(field, Collections.singletonList(coord));
	}

	public AffinePoint(Field<T> field, T coord1, T coord2) {
		this.coords = new ArrayList<T>();
		this.coords.add(coord1);
		this.coords.add(coord2);
		this.init(field);
	}

	public AffinePoint(Field<T> field, T coord1, T coord2, T coord3) {
		this.coords = new ArrayList<T>();
		this.coords.add(coord1);
		this.coords.add(coord2);
		this.coords.add(coord3);
		this.init(field);
	}

	public AffinePoint(Field<T> field, List<T> coords) {
		this.coords = new ArrayList<T>();
		this.coords.addAll(coords);
		this.init(field);
	}

	private void init(Field<T> field) {
		this.dim = this.coords.size();
	}

	public int getDim() {
		return dim;
	}

	public List<T> getCoords() {
		return Collections.unmodifiableList(coords);
	}

	public T getCoord(int i) {
		return this.coords.get(i - 1);
	}

	public Ideal<Polynomial<T>> asIdeal(PolynomialRing<T> ring) {
		List<Polynomial<T>> generators = new ArrayList<>();
		for (int i = 0; i < coords.size(); i++) {
			T t = coords.get(i);
			generators.add(ring.subtract(ring.getVar(i + 1), ring.getEmbedding(t)));
		}
		return ring.getIdeal(generators);
	}

	@Override
	public boolean equals(Object O) {
		if (!(O instanceof AffinePoint))
			return false;
		@SuppressWarnings("unchecked")
		AffinePoint<T> T = (AffinePoint<T>) O;
		if (T.getDim() != this.getDim())
			return false;
		for (int i = 1; i <= this.getDim(); i++) {
			if (!this.getCoord(i).equals(T.getCoord(i)))
				return false;
		}
		return true;
	}

	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append("(");
		boolean first = true;
		for (T c : this.coords) {
			if (first)
				first = false;
			else
				buf.append(", ");
			buf.append(c);
		}
		buf.append(")");
		return buf.toString();
	}

	@Override
	public int compareTo(AffinePoint<T> p) {
		for (int i = 1; i <= this.getDim(); i++) {
			int cmp = this.getCoord(i).compareTo(p.getCoord(i));
			if (cmp != 0)
				return cmp;
		}
		return 0;
	}

}
