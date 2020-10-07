package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.polynomials.Monomial;
import fields.vectors.Vector;

public class CoordinateRing<T extends Element<T>> extends AbstractAlgebra<Polynomial<T>, CoordinateRingElement<T>>
		implements Ring<CoordinateRingElement<T>> {
	private PolynomialRing<T> ring;
	private Ideal<Polynomial<T>> ideal;
	private int dimension;

	public static class CoordinateRingElement<T extends Element<T>> extends AbstractElement<CoordinateRingElement<T>> {
		private Polynomial<T> polynomial;

		public CoordinateRingElement(Ideal<Polynomial<T>> ideal, Polynomial<T> polynomial) {
			this.polynomial = ideal.residue(polynomial);
		}

		private CoordinateRingElement(Ideal<Polynomial<T>> ideal, Polynomial<T> polynomial, boolean reduce) {
			if (reduce) {
				this.polynomial = ideal.residue(polynomial);
			} else {
				this.polynomial = polynomial;
			}
		}

		public String toString() {
			return this.polynomial.toString();
		}

		public Polynomial<T> getElement() {
			return this.polynomial;
		}

		@Override
		public int compareTo(CoordinateRingElement<T> o) {
			return this.polynomial.compareTo(o.polynomial);
		}
	}

	public CoordinateRing(PolynomialRing<T> ring, Ideal<Polynomial<T>> ideal) {
		this.ring = ring;
		this.ideal = ideal;
		int i = 0;
		int[][] leadingMonomials = new int[this.ideal.generators().size()][this.ring.numberOfVariables()];
		for (Polynomial<T> generator : this.ideal.generators()) {
			leadingMonomials[i] = generator.leadingMonomial().exponents();
			i++;
		}

		this.dimension = this.ring.numberOfVariables()
				- this.makeset(0, new int[this.ring.numberOfVariables() + 1], leadingMonomials)[this.ring.numberOfVariables()];
	}

	private int[] makeset(int level, int[] set, int[][] leadingMonomials) {
		if (level == this.ideal.generators().size())
			return set;
		int[] optimalset = new int[set.length];
		optimalset[set.length - 1] = -1;
		for (int i = 0; i < leadingMonomials[level].length; i++) {
			if (leadingMonomials[level][i] == 0)
				continue;
			if (set[i] == 1)
				return this.makeset(level + 1, set, leadingMonomials);
			set[i] = 1;
			set[set.length - 1]++;
			int[] sethere = this.makeset(level + 1, set, leadingMonomials);
			if (optimalset[set.length - 1] == -1 || optimalset[set.length - 1] > sethere[set.length - 1])
				optimalset = Arrays.copyOf(sethere, sethere.length);
			set[i] = 0;
			set[set.length - 1]--;
		}
		return optimalset;
	}

	@Override
	public String toString() {
		return this.ring.toString() + "/" + this.ideal.toString();
	}

	public PolynomialRing<T> getPolynomialRing() {
		return ring;
	}

	public Ring<Polynomial<T>> getRing() {
		return ring;
	}

	public Ideal<Polynomial<T>> getIdeal() {
		return ideal;
	}

	@Override
	public int krullDimension() {
		return dimension;
	}

	public boolean isFree() {
		throw new ArithmeticException();
	}

	public List<CoordinateRingElement<T>> getGenerators() {
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		for (int i = 1; i <= this.ring.numberOfVariables(); i++) {
			result.add(getVar(i));
		}
		return result;
	}

	public CoordinateRingElement<T> getVar(int var) {
		return this.getEmbedding(this.ring.getVar(var));
	}

	public CoordinateRingElement<T> getEmbedding(Polynomial<T> t) {
		return new CoordinateRingElement<>(this.ideal, t);
	}

	private CoordinateRingElement<T> getEmbedding(Polynomial<T> t, boolean reduce) {
		return new CoordinateRingElement<>(this.ideal, t, reduce);
	}

	public CoordinateRingElement<T> getEmbedding(T t) {
		return this.getEmbedding(this.ring.getEmbedding(t), false);
	}

	@Override
	public CoordinateRingElement<T> zero() {
		return this.getEmbedding(this.ring.zero(), false);
	}

	@Override
	public CoordinateRingElement<T> one() {
		return this.getEmbedding(this.ring.one(), false);
	}

	@Override
	public CoordinateRingElement<T> add(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.add(t1.polynomial, t2.polynomial), false);
	}

	@Override
	public CoordinateRingElement<T> negative(CoordinateRingElement<T> t) {
		return this.getEmbedding(this.ring.negative(t.polynomial), false);
	}

	@Override
	public CoordinateRingElement<T> multiply(CoordinateRingElement<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1.polynomial, t2.polynomial));
	}

	public CoordinateRingElement<T> multiply(T t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> scalarMultiply(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	@Override
	public CoordinateRingElement<T> inverse(CoordinateRingElement<T> t) {
		PolynomialRing<T> r = ring.addVariableWithElimination(1);
		Polynomial<T> s = r.getVar(1);
		Polynomial<T> f = r.getEmbeddingWithElimination(t.getElement(), 1);
		Polynomial<T> inv = r.subtract(r.one(), r.multiply(s, f));
		List<Polynomial<T>> basis = new ArrayList<Polynomial<T>>();
		basis.add(inv);
		for (Polynomial<T> b : ideal.generators()) {
			basis.add(r.getEmbeddingWithElimination(b, 1));
		}
		Ideal<Polynomial<T>> i = r.getIdeal(basis);
		Polynomial<T> inverted = i.generators().get(0);
		int[] expectedExponents = new int[r.numberOfVariables()];
		expectedExponents[0] = 1;
		if (!inverted.leadingMonomial().equals(r.getMonomial(expectedExponents))) {
			throw new ArithmeticException("No inverse found!");
		}
		Polynomial<T> shifted = ring.getEmbeddingWithElimination(r.add(s, r.negative(inverted)), -1);
		return this.getEmbedding(shifted);
	}

	@Override
	public boolean isUnit(CoordinateRingElement<T> t) {
		Ideal<Polynomial<T>> tIdeal = this.ring.getIdeal(Collections.singletonList(t.polynomial));
		Ideal<Polynomial<T>> degenerate = this.ring.getIdeal(Collections.singletonList(this.ring.one()));
		return this.ring.add(this.ideal, tIdeal).equals(degenerate);
	}
	
	public CoordinateRingElement<T> projectToUnit(CoordinateRingElement<T> t) {
		return getEmbedding(ring.projectToUnit(t.polynomial));
	}

	@Override
	public CoordinateRingElement<T> getRandomElement() {
		throw new UnsupportedOperationException();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<CoordinateRingElement<T>> iterator() {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isFinite() {
		return this.ring.getRing().isFinite() && this.krullDimension() == 0;
	}

	@Override
	public BigInteger characteristic() {
		return this.ring.characteristic();
	}

	@Override
	public boolean isIntegral() {
		return false; // TODO check whether integral
	}

	@Override
	public boolean isZeroDivisor(CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isDivisible(CoordinateRingElement<T> dividend, CoordinateRingElement<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> quotientAndRemainder(CoordinateRingElement<T> dividend,
			CoordinateRingElement<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public BigInteger euclidMeasure(CoordinateRingElement<T> t) {
		return null;
	}

	@Override
	public Iterable<CoordinateRingElement<T>> getUnits() {
		throw new UnsupportedOperationException();
	}

	@Override
	public Ideal<CoordinateRingElement<T>> getIdeal(List<CoordinateRingElement<T>> generators) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Ideal<CoordinateRingElement<T>> intersect(Ideal<CoordinateRingElement<T>> t1,
			Ideal<CoordinateRingElement<T>> t2) {
		throw new UnsupportedOperationException();
	}

	public CoordinateRingElement<T> substitute(CoordinateRingElement<T> t, List<CoordinateRingElement<T>> values) {
		CoordinateRingElement<T> result = this.zero();
		for (Monomial m : t.polynomial.monomials()) {
			CoordinateRingElement<T> value = this.getEmbedding(t.polynomial.coefficient(m));
			int exponents[] = m.exponents();
			for (int i = 0; i < exponents.length; i++) {
				value = this.multiply(value, this.power(values.get(i), exponents[i]));
			}
			result = this.add(result, value);
		}
		return result;
	}

	@Override
	public boolean isGeneratingAlgebra(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> getAlgebraGenerators() {
		return Collections.singletonList(one());
	}

	@Override
	public boolean isLinearIndependent(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isGeneratingModule(List<CoordinateRingElement<T>> s) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<CoordinateRingElement<T>> getModuleGenerators() {
		return Collections.singletonList(one());
	}

	@Override
	public Vector<Polynomial<T>> asVector(CoordinateRingElement<T> s) {
		return new Vector<>(Collections.singletonList(s.polynomial));
	}
	
	public static class CoordinateIdeal<T extends Element<T>> extends AbstractIdeal<CoordinateRingElement<T>> {
		private CoordinateRing<T> ring;
		private List<CoordinateRingElement<T>> generators;
		private Ideal<Polynomial<T>> asPolynomialIdeal;
		
		private CoordinateIdeal(CoordinateRing<T> ring, List<CoordinateRingElement<T>> generators) {
			super(ring);
			List<Polynomial<T>> generatorPolynomials = new ArrayList<>();
			for (CoordinateRingElement<T> generator : generators) {
				generatorPolynomials.add(generator.polynomial);
			}
			generatorPolynomials.addAll(ring.getIdeal().generators());
			this.ring = ring;
			this.asPolynomialIdeal = ring.getPolynomialRing().getIdeal(generatorPolynomials);
			this.generators = new ArrayList<>();
			for (Polynomial<T> generator : this.asPolynomialIdeal.generators()) {
				Polynomial<T> reduced = ring.getIdeal().residue(generator);
				if (!reduced.equals(ring.getPolynomialRing().zero())) {
					this.generators.add(ring.getEmbedding(reduced));
				}
			}
		}
		
		@Override
		public boolean isPrime() {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public boolean isMaximal() {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public List<CoordinateRingElement<T>> generators() {
			return Collections.unmodifiableList(generators);
		}

		@Override
		public List<CoordinateRingElement<T>> generate(CoordinateRingElement<T> t) {
			List<Polynomial<T>> generate = this.asPolynomialIdeal.generate(t.polynomial);
			List<Polynomial<T>> polynomialGenerators = this.asPolynomialIdeal.generators();
			List<CoordinateRingElement<T>> result = new ArrayList<>();
			for (int i = 0; i < polynomialGenerators.size(); i++) {
				Polynomial<T> generator = polynomialGenerators.get(i);
				Polynomial<T> reduced = ring.getIdeal().residue(generator);
				if (!reduced.equals(ring.getPolynomialRing().zero())) {
					result.add(ring.getEmbedding(generate.get(i)));
				}
			}
			return result;
		}

		@Override
		public CoordinateRingElement<T> residue(CoordinateRingElement<T> t) {
			return ring.getEmbedding(asPolynomialIdeal.residue(t.polynomial));
		}

		@Override
		public boolean contains(CoordinateRingElement<T> t) {
			return asPolynomialIdeal.contains(t.polynomial);
		}

		@Override
		public boolean isFinite() {
			return ring.isFinite();
		}

		@Override
		public BigInteger getNumberOfElements() throws InfinityException {
			throw new UnsupportedOperationException();
		}
		
	}
}
