package varieties;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import fields.Element;
import fields.Field;
import fields.Polynomial;
import fields.PolynomialRing;

public class ProjectivePoint<T extends Element> implements Element {
	private Field<T> field;
	private int dim;
	private List<T> coords;
	private int nonzero;
	private PolynomialRing<T>.Ideal ideal;

        @SafeVarargs
	public ProjectivePoint(Field<T> field, T... coords) {
		this.coords = new ArrayList<T>();
		this.coords.addAll(Arrays.<T>asList(coords));
		this.init(field);
	}
	public ProjectivePoint(Field<T> field, List<T> coords) {
		this.coords = new ArrayList<T>();
		this.coords.addAll(coords);
		this.init(field);
	}
	private void init(Field<T> field) {
		this.field = field;
		boolean nonzerofound = false;
		int i = 1;
		for (T coord : this.coords) {
			if (!coord.equals(this.field.zero())) {
				nonzerofound = true;
				nonzero = i;
			}
			i++;
		}
		if (!nonzerofound)
			throw new ArithmeticException("Divison of zero by zero!");
		this.dim = this.coords.size() - 1;
	}
	public PolynomialRing<T>.Ideal asIdeal(PolynomialRing<T> ring) {
		if (this.ideal != null)
			return this.ideal;
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
		int nonzero = this.getNonZero() - 1;
		List<T> coeff = new ArrayList<T>();
		for (int i = 0; i < this.dim + 1; i++)
			coeff.add(this.field.zero());
		for (int i = 0; i < this.dim + 1; i++) {
			if (i == nonzero)
				continue;
			coeff.set(i, this.field.one());
			coeff.set(nonzero, this.field.negative(this.getDehomogenisedCoord(i + 1, nonzero + 1)));
			list.add(ring.getLinear(coeff));
			coeff.set(i, this.field.zero());
		}
		this.ideal = ring.getIdeal(list);
		return this.asIdeal(ring);
	}
	public int getDim() {
		return dim;
	}
	public int getNonZero() {
		return this.nonzero;
	}
	public List<T> getCoords() {
		return Collections.unmodifiableList(coords);
	}
	public T getCoord(int i) {
		return this.coords.get(i - 1);
	}
	public T getDehomogenisedCoord(int i, int dehomogencoord) {
		if (this.getCoord(dehomogencoord).equals(this.field.zero()))
			throw new ArithmeticException("Divison by zero");
		else
			return field.divide(this.getCoord(i), this.getCoord(dehomogencoord));
	}
	public AffinePoint<T> getDehomogenous(int dehomogencoord) {
		List<T> coord = new ArrayList<T>();
		for (int i = 1; i <= this.getDim() + 1; i++)
			coord.add(this.getDehomogenisedCoord(i, dehomogencoord));
		return new AffinePoint<T>(field, coord);
	}
	@Override
	public boolean equals(Object O) {
		if (!(O instanceof ProjectivePoint))
			return false;
		@SuppressWarnings("unchecked")
		ProjectivePoint<T> T = (ProjectivePoint<T>)O;
		if (T.getDim() != this.getDim())
			return false;
		if (this.field.zero().equals(T.getCoord(nonzero)))
			return false;
		for (int i = 1; i <= this.getDim() + 1; i++) {
			if (i == this.nonzero)
				continue;
			if (!this.getDehomogenisedCoord(i, this.nonzero).equals(T.getDehomogenisedCoord(i, this.nonzero)))
				return false;
		}
		return true;
	}
	@Override
	public int hashCode() {
		int res = 37;
		for (int i = 1; i <= this.getDim() + 1; i++) {
			if (i == this.nonzero)
				continue;
			res += 17 * this.getDehomogenisedCoord(i, this.nonzero).hashCode();
		}
		return res;
	}
	@Override
	public String toString() {
		StringBuffer buf = new StringBuffer();
		buf.append("[");
		boolean first = true;
		for (int i = 1; i <= this.dim + 1; i++) {
			if (first)
				first = false;
			else
				buf.append(":");
			buf.append(this.getDehomogenisedCoord(i, nonzero));
		}
		buf.append("]");
		return buf.toString();
	}
	@Override
	public int compareTo(Element o) {
		@SuppressWarnings("unchecked")
		ProjectivePoint<T> p = (ProjectivePoint<T>)o;
		if (p.getDim() != this.getDim())
			return this.getDim() - p.getDim();
		if (this.nonzero != p.nonzero)
			return this.nonzero - p.nonzero;
		for (int i = 1; i <= this.getDim() + 1; i++) {
			if (i == this.nonzero)
				continue;
			int cmp = this.getDehomogenisedCoord(i, this.nonzero).compareTo(p.getDehomogenisedCoord(i, this.nonzero)); 
			if (cmp != 0)
				return cmp;
		}
		return 0;
	}
}
