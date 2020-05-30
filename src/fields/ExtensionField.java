package fields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.ExtensionField.ExtensionFieldElement;
import fields.Polynomial.Monomial;

public class ExtensionField<T extends Element> extends AbstractField<ExtensionFieldElement<T>>
		implements Field<ExtensionFieldElement<T>> {
	private Polynomial<T> minimalpolynomial;
	private Field<T> field;
	private int degree;
	private PolynomialRing<T> ring;
	private Map<T, ExtensionFieldElement<T>> embed;
	private Comparator<Polynomial.Monomial> comparator;

	public ExtensionField(Polynomial<T> minimalpolynomial, Field<T> field) {
		this.comparator = Polynomial.LEX;
		this.field = field;
		this.ring = minimalpolynomial.getRing();
		this.minimalpolynomial = minimalpolynomial.normalize();
		if (this.minimalpolynomial.getNumVars() != 1)
			throw new ArithmeticException("Not univariate!");
		if (this.ring.factorization(this.minimalpolynomial).size() != 1) {
			throw new ArithmeticException("Not irreducible!");
		}
		this.degree = this.minimalpolynomial.getDegree();
		this.embed = new TreeMap<T, ExtensionFieldElement<T>>();
	}

	private static <T extends Element> Polynomial<T> findIrreduciblePolynomial(int degree, Field<T> field) {
		Polynomial<T> result;
		PolynomialRing<T> r = new PolynomialRing<T>(field, 1, Polynomial.LEX);
		do {
			Map<Polynomial.Monomial, T> coeffs = new TreeMap<>(r.getComparator());
			coeffs.put(new Polynomial.Monomial(new int[] { degree }), field.one());
			for (int i = 0; i < degree; i++) {
				coeffs.put(new Polynomial.Monomial(new int[] { i }), field.getRandomElement());
			}
			result = new Polynomial<T>(coeffs, r);
		} while (r.factorization(result).size() != 1);
		return result;
	}

	public ExtensionField(int degree, Field<T> field) {
		this(findIrreduciblePolynomial(degree, field), field);
	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	@Override
	public String toString() {
		return this.field.toString() + "[X]/(" + this.minimalpolynomial.toString() + ")";
	}

	public ExtensionFieldElement<T> getEmbedding(T t) {
		if (this.embed.containsKey(t))
			return this.embed.get(t);
		List<T> list = new ArrayList<T>();
		list.add(t);
		for (int i = 1; i < this.degree; i++)
			list.add(this.field.zero());
		ExtensionFieldElement<T> e = new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
		this.embed.put(t, e);
		return this.getEmbedding(t);
	}

	@Override
	public ExtensionFieldElement<T> zero() {
		return this.getEmbedding(this.field.zero());
	}

	@Override
	public ExtensionFieldElement<T> one() {
		return this.getEmbedding(this.field.one());
	}
	
	public ExtensionFieldElement<T> alpha() {
		return fromPolynomial(ring.getVar(1));
	}

	@Override
	public BigInteger characteristic() {
		return this.field.characteristic();
	}

	@Override
	public ExtensionFieldElement<T> add(ExtensionFieldElement<T> t1, ExtensionFieldElement<T> t2) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.add(t1.getElement().get(i), t2.getElement().get(i)));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}

	@Override
	public ExtensionFieldElement<T> negative(ExtensionFieldElement<T> t) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.negative(t.getElement().get(i)));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}

	@Override
	public ExtensionFieldElement<T> multiply(ExtensionFieldElement<T> t1, ExtensionFieldElement<T> t2) {
		Polynomial<T> res = this.ring.multiply(this.asPolynomial(t1), this.asPolynomial(t2));
		return this.fromPolynomial(this.ring.quotientAndRemainder(res, this.minimalpolynomial).get(1));
	}

	@Override
	public ExtensionFieldElement<T> inverse(ExtensionFieldElement<T> t) {
		List<Polynomial<T>> egcd = this.ring.extendedEuclidean(this.asPolynomial(t), this.minimalpolynomial);
		if (!this.ring.isUnit(egcd.get(0))) {
			throw new ArithmeticException("Division by zero");
		}
		ExtensionFieldElement<T> inv = this
				.fromPolynomial(this.ring.multiply(egcd.get(1), this.ring.inverse(egcd.get(0))));
		return inv;
	}

	public Polynomial<T> asPolynomial(ExtensionFieldElement<T> t) {
		Map<Monomial, T> map = new TreeMap<Monomial, T>(this.comparator);
		for (int i = 0; i < this.degree; i++)
			map.put(new Monomial(new int[] { i }), t.getElement().get(i));
		return new Polynomial<T>(map, this.ring);
	}

	public Vector<T> asVector(ExtensionFieldElement<T> t) {
		return new Vector<T>(t.coeff);
	}
	
	public ExtensionFieldElement<T> fromVector(Vector<T> t) {
		return new ExtensionFieldElement<>(field, minimalpolynomial, t.asList());
	}

	public Matrix<T> asMatrix(ExtensionFieldElement<T> t) {
		List<List<T>> m = new ArrayList<>();
		Vector<T> asVector = asVector(t);
		List<T> minPolyList = new ArrayList<>();
		for (int i = 0; i < this.degree; i++) {
			m.add(new ArrayList<>());
			minPolyList.add(minimalpolynomial.getCoefficient(new Polynomial.Monomial(new int[] { i })));
		}
		Vector<T> minPoly = new Vector<T>(minPolyList);
		FreeModule<T> space = new FreeModule<>(field, degree);
		for (int j = 0; j < this.degree; j++) {
			List<T> nextVector = new ArrayList<>();
			nextVector.add(field.zero());
			for (int i = 0; i < this.degree; i++) {
				m.get(i).add(asVector.get(i + 1));
				nextVector.add(asVector.get(i + 1));
			}
			T high = nextVector.get(degree);
			asVector = space.subtract(new Vector<>(nextVector.subList(0, degree)), space.scalarMultiply(high, minPoly));
		}
		return new Matrix<T>(m);
	}

	private ExtensionFieldElement<T> fromPolynomial(Polynomial<T> t) {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(t.getCoefficient(new Monomial(new int[] { i })));
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}

	@Override
	public ExtensionFieldElement<T> getRandomElement() {
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(this.field.getRandomElement());
		return new ExtensionFieldElement<T>(this.field, this.minimalpolynomial, list);
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		BigInteger num = this.field.getNumberOfElements();
		return num.pow(this.degree);
	}

	@Override
	public Iterator<ExtensionFieldElement<T>> iterator() {
		return new Iterator<ExtensionFieldElement<T>>() {
			private List<Iterator<T>> it = null;
			private List<T> list = new ArrayList<T>();

			private void init() {
				if (this.it != null)
					return;
				this.it = new ArrayList<Iterator<T>>();
				for (int i = 0; i < degree; i++) {
					this.it.add(field.iterator());
					if (i == 0)
						this.list.add(null);
					else
						this.list.add(this.it.get(i).next());
				}
			}

			@Override
			public boolean hasNext() {
				init();
				for (int i = 0; i < degree; i++) {
					if (this.it.get(i).hasNext())
						return true;
				}
				return false;
			}

			@Override
			public ExtensionFieldElement<T> next() {
				init();
				boolean broken = false;
				for (int i = 0; i < degree; i++) {
					if (this.it.get(i).hasNext()) {
						this.list.set(i, this.it.get(i).next());
						for (int j = 0; j < i; j++) {
							this.it.set(j, field.iterator());
							this.list.set(j, this.it.get(j).next());
						}
						broken = true;
						break;
					}
				}
				if (!broken)
					throw new RuntimeException();
				return new ExtensionFieldElement<T>(field, minimalpolynomial, list);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}

		};
	}

	public static class ExtensionFieldElement<T extends Element> implements Element {
		private Field<T> base;
		private List<T> coeff;

		public ExtensionFieldElement(Field<T> base, Polynomial<T> minimalpolynomial, List<T> coeff) {
			this.base = base;
			this.coeff = new ArrayList<T>();
			if (minimalpolynomial.getNumVars() != 1 || coeff.size() != minimalpolynomial.getDegree())
				throw new ArithmeticException("Invalid minimal polynomial");
			this.coeff.addAll(coeff);
		}

		public boolean equals(Object o) {
			if (!(o instanceof ExtensionField.ExtensionFieldElement))
				return false;
			@SuppressWarnings("unchecked")
			ExtensionFieldElement<T> t = (ExtensionFieldElement<T>) o;
			return this.coeff.equals(t.coeff);
		}

		public String toString() {
			boolean first = true;
			boolean onlyzero = true;
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < this.coeff.size(); i++) {
				if (!this.coeff.get(i).equals(this.base.zero())) {
					onlyzero = false;
					if (first)
						first = false;
					else
						buf.append(" + ");
					if (!this.coeff.get(i).equals(this.base.one()) || i == 0)
						buf.append(this.coeff.get(i));
					if (i == 1)
						buf.append("a");
					if (i > 1)
						buf.append("a^" + i);
				}
			}
			if (onlyzero)
				return "0";
			return buf.toString();
		}

		public List<T> getElement() {
			return Collections.unmodifiableList(this.coeff);
		}

		@Override
		public int compareTo(Element e) {
			@SuppressWarnings("unchecked")
			ExtensionFieldElement<T> o = (ExtensionFieldElement<T>) e;
			for (int i = 0; i < this.coeff.size(); i++) {
				int cmp = this.coeff.get(i).compareTo(o.coeff.get(i));
				if (cmp != 0)
					return cmp;
			}
			return 0;
		}
	}

	public static <T extends Element> PolynomialRing<T> getPolynomialRing(Field<T> base) {
		return new PolynomialRing<T>(base, 1, Polynomial.LEX);
	}
}
