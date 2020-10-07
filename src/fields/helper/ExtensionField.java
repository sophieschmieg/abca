package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.interfaces.Algebra;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.VectorSpace;
import fields.polynomials.Monomial;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class ExtensionField<T extends Element<T>> extends AbstractField<ExtensionFieldElement<T>>
		implements VectorSpace<T, ExtensionFieldElement<T>>, Algebra<T, ExtensionFieldElement<T>> {
	private Polynomial<T> minimalpolynomial;
	private Field<T> field;
	private int degree;
	private PolynomialRing<T> ring;
	private FreeModule<T> asVectorSpace;
	private MatrixAlgebra<T> algebra;
	private Map<T, ExtensionFieldElement<T>> embed;
	
	public ExtensionField(Polynomial<T> minimalpolynomial, Field<T> field) {
		this(minimalpolynomial, field, false);
	}

	public ExtensionField(Polynomial<T> minimalpolynomial, Field<T> field, boolean checkIrreducible) {
		this.field = field;
		this.ring = field.getUnivariatePolynomialRing();
		this.minimalpolynomial = ring.normalize(minimalpolynomial);
		if (checkIrreducible && field.factorization(this.minimalpolynomial).size() != 1) {
			throw new ArithmeticException("Not irreducible!");
		}
		this.degree = this.minimalpolynomial.degree();
		this.embed = new TreeMap<T, ExtensionFieldElement<T>>();
		this.asVectorSpace = new FreeModule<>(this.field, this.degree);
		this.algebra = asVectorSpace.matrixAlgebra();
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

	public Polynomial<T> minimalPolynomial() {
		return minimalpolynomial;
	}

	public int degree() {
		return this.degree;
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
		Map<Monomial, T> map = new TreeMap<Monomial, T>();
		for (int i = 0; i < this.degree; i++)
			map.put(ring.getMonomial(new int[] { i }), t.getElement().get(i));
		return ring.getPolynomial(map);
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
			minPolyList.add(minimalpolynomial.coefficient(ring.getMonomial(new int[] { i })));
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

	public T norm(ExtensionFieldElement<T> t) {
		return algebra.determinant(asMatrix(t));
	}

	public T trace(ExtensionFieldElement<T> t) {
		return algebra.trace(asMatrix(t));
	}

	public Polynomial<T> minimalPolynomial(ExtensionFieldElement<T> t) {
		Matrix<T> matrixElement = asMatrix(t);
		List<List<Polynomial<T>>> matrix = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			List<Polynomial<T>> row = new ArrayList<>();
			matrix.add(row);
			for (int j = 0; j < degree; j++) {
				Polynomial<T> entry = ring.getEmbedding(matrixElement.entry(i + 1, j + 1));
				if (i == j) {
					entry = ring.subtract(entry, ring.getVar(1));
				}
				row.add(entry);
			}
		}
		Matrix<Polynomial<T>> m = new Matrix<>(matrix);
		MatrixAlgebra<Polynomial<T>> algebra = new FreeModule<>(ring, degree).matrixAlgebra();
		return ring.normalize(algebra.determinant(m));
	}

	public ExtensionFieldElement<T> fromPolynomial(Polynomial<T> t) {
		if (t.degree() >= this.degree) {
			t = ring.quotientAndRemainder(t, minimalpolynomial).get(1);
		}
		List<T> list = new ArrayList<T>();
		for (int i = 0; i < this.degree; i++)
			list.add(t.coefficient(ring.getMonomial(new int[] { i })));
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

	public static class ExtensionFieldElement<T extends Element<T>> implements Element<ExtensionFieldElement<T>> {
		private Field<T> base;
		private List<T> coeff;

		public ExtensionFieldElement(Field<T> base, Polynomial<T> minimalpolynomial, List<T> coeff) {
			this.base = base;
			this.coeff = new ArrayList<T>();
			if (minimalpolynomial.numberOfVariables() != 1 || coeff.size() != minimalpolynomial.degree())
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
		public int compareTo(ExtensionFieldElement<T> o) {
			for (int i = 0; i < this.coeff.size(); i++) {
				int cmp = this.coeff.get(i).compareTo(o.coeff.get(i));
				if (cmp != 0)
					return cmp;
			}
			return 0;
		}
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(T t, ExtensionFieldElement<T> s) {
		return multiply(getEmbedding(t), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(int n, ExtensionFieldElement<T> s) {
		return multiply(getInteger(n), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(BigInteger n, ExtensionFieldElement<T> s) {
		return multiply(getInteger(n), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(T t1, T t2, ExtensionFieldElement<T> s) {
		return multiply(getEmbedding(t1), getEmbedding(t2), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(int n, T t, ExtensionFieldElement<T> s) {
		return multiply(getInteger(n), getEmbedding(t), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(BigInteger n, T t, ExtensionFieldElement<T> s) {
		return multiply(getInteger(n), getEmbedding(t), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(int n, T t1, T t2, ExtensionFieldElement<T> s) {
		return multiply(multiply(getInteger(n), getEmbedding(t1)), getEmbedding(t2), s);
	}

	@Override
	public ExtensionFieldElement<T> scalarMultiply(BigInteger n, T t1, T t2, ExtensionFieldElement<T> s) {
		return multiply(multiply(getInteger(n), getEmbedding(t1)), getEmbedding(t2), s);
	}

	public boolean isGeneratingAlgebra(List<ExtensionFieldElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Ring<T> getRing() {
		return field;
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public boolean isLinearIndependent(List<ExtensionFieldElement<T>> s) {
	  List<Vector<T>> asVectors = new ArrayList<>();
		for (ExtensionFieldElement<T> element : s) {
		  asVectors.add(asVector(element));
	  }
		return asVectorSpace.isLinearIndependent(asVectors);
	}

	@Override
	public boolean isGeneratingModule(List<ExtensionFieldElement<T>> s) {
		  List<Vector<T>> asVectors = new ArrayList<>();
			for (ExtensionFieldElement<T> element : s) {
			  asVectors.add(asVector(element));
		  }
			return asVectorSpace.isGeneratingModule(asVectors);
		}

	@Override
	public List<ExtensionFieldElement<T>> getModuleGenerators() {
		List<ExtensionFieldElement<T>> generators = new ArrayList<>();
		ExtensionFieldElement<T> generator = one();
		for (int i = 0; i < degree; i++) {
			generators.add(generator);
			generator = multiply(generator, alpha());
		}
		return generators;
	}

	@Override
	public List<ExtensionFieldElement<T>> getAlgebraGenerators() {
		return Collections.singletonList(alpha());
	}

	@Override
	public Field<T> getField() {
		return field;
	}

	@Override
	public List<ExtensionFieldElement<T>> getBasis() {
		return getModuleGenerators();
	}

}
