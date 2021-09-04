package fields.helper;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Vector;

public class CoordinateRing<T extends Element<T>> extends AbstractAlgebra<Polynomial<T>, CoordinateRingElement<T>>
		implements Ring<CoordinateRingElement<T>> {
	private PolynomialRing<T> ring;
	private PolynomialIdeal<T> ideal;
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

	public static <T extends Element<T>> CoordinateRing<T> tensorProduct(CoordinateRing<T> t1, CoordinateRing<T> t2) {
		return tensorProduct(t1, t2, Monomial.GREVLEX);
	}

	public static <T extends Element<T>> CoordinateRing<T> tensorProduct(CoordinateRing<T> t1, CoordinateRing<T> t2,
			Comparator<Monomial> comparator) {
		if (!t1.getPolynomialRing().getRing().equals(t2.getPolynomialRing().getRing())) {
			throw new ArithmeticException("Not over common base!");
		}
		Ring<T> base = t1.getPolynomialRing().getRing();
		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(base,
				t1.getPolynomialRing().numberOfVariables() + t2.getPolynomialRing().numberOfVariables(), comparator);
		int[] map1 = new int[t1.getPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map1.length; i++) {
			map1[i] = i;
		}
		int[] map2 = new int[t2.getPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map2.length; i++) {
			map2[i] = i + map1.length;
		}
		List<Polynomial<T>> polynomials = new ArrayList<>();
		for (Polynomial<T> polynomial : t1.getIdeal().generators()) {
			polynomials.add(polynomialRing.getEmbedding(polynomial, map1));
		}
		for (Polynomial<T> polynomial : t2.getIdeal().generators()) {
			polynomials.add(polynomialRing.getEmbedding(polynomial, map2));
		}
		return new CoordinateRing<>(polynomialRing, polynomialRing.getIdeal(polynomials));
	}

	public CoordinateRing(CoordinateRing<T> ring, CoordinateIdeal<T> ideal) {
		this(ring.getPolynomialRing(), ideal.asPolynomialIdeal());
	}

	public CoordinateRing(PolynomialRing<T> ring, PolynomialIdeal<T> ideal) {
		this.ring = ring;
		this.ideal = ideal;
		int i = 0;
		int[][] leadingMonomials = new int[this.ideal.generators().size()][this.ring.numberOfVariables()];
		for (Polynomial<T> generator : this.ideal.generators()) {
			leadingMonomials[i] = generator.leadingMonomial().exponents();
			i++;
		}
		if (ideal.contains(ring.one())) {
			this.dimension = -1;
		} else {
			this.dimension = this.ring.numberOfVariables() - this.makeset(0, new int[this.ring.numberOfVariables() + 1],
					leadingMonomials)[this.ring.numberOfVariables()];
		}
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
	public Exactness exactness() {
		return ring.exactness();
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

	public PolynomialIdeal<T> getIdeal() {
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

	@Override
	public CoordinateRingElement<T> getEmbedding(Polynomial<T> t) {
		return new CoordinateRingElement<>(this.ideal, ring.getEmbedding(t));
	}

	private CoordinateRingElement<T> getEmbedding(Polynomial<T> t, boolean reduce) {
		return new CoordinateRingElement<>(this.ideal, ring.getEmbedding(t), reduce);
	}

	public CoordinateRingElement<T> getRingEmbedding(T t) {
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
		return this.getEmbedding(this.ring.add(t1.polynomial, t2.polynomial));
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

	public CoordinateRingElement<T> multiplyAssign(T t1, CoordinateRingElement<T> t2) {
		return getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> scalarMultiply(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return this.getEmbedding(this.ring.multiply(t1, t2.polynomial));
	}

	public CoordinateRingElement<T> scalarMultiplyAssign(Polynomial<T> t1, CoordinateRingElement<T> t2) {
		return getEmbedding(this.ring.multiply(t1, t2.polynomial));
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
		for (Polynomial<T> generator : ideal.generators()) {
			if (!ring.isIrreducible(generator)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isZeroDivisor(CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isCommutative() {
		return ring.isCommutative();
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<CoordinateRingElement<T>> uniqueFactorization(CoordinateRingElement<T> t) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(CoordinateRingElement<T> dividend, CoordinateRingElement<T> divisor) {
		throw new UnsupportedOperationException();
	}

	@Override
	public QuotientAndRemainderResult<CoordinateRingElement<T>> quotientAndRemainder(CoordinateRingElement<T> dividend,
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
	public IdealResult<CoordinateRingElement<T>, CoordinateIdeal<T>> getIdealWithTransforms(
			List<CoordinateRingElement<T>> generators) {
		throw new UnsupportedOperationException();
	}

	public CoordinateIdeal<T> getIdeal(List<CoordinateRingElement<T>> generators) {
		return new CoordinateIdeal<>(this, generators);
	}

	public CoordinateIdeal<T> getIdeal(Ideal<Polynomial<T>> ideal) {
		List<CoordinateRingElement<T>> generators = new ArrayList<>();
		for (Polynomial<T> t : ideal.generators()) {
			generators.add(getEmbedding(t));
		}
		return getIdeal(generators);
	}

	@Override
	public Ideal<CoordinateRingElement<T>> intersect(Ideal<CoordinateRingElement<T>> t1,
			Ideal<CoordinateRingElement<T>> t2) {
		throw new UnsupportedOperationException();
	}

	@Override
	public Ideal<CoordinateRingElement<T>> radical(Ideal<CoordinateRingElement<T>> t) {
		throw new UnsupportedOperationException();
	}

	public CoordinateRingElement<T> substitute(CoordinateRingElement<T> t, List<CoordinateRingElement<T>> values) {
		CoordinateRingElement<T> result = this.zero();
		for (Monomial m : t.polynomial.monomials()) {
			CoordinateRingElement<T> value = this.getRingEmbedding(t.polynomial.coefficient(m));
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
	public List<List<Polynomial<T>>> nonTrivialCombinations(List<CoordinateRingElement<T>> s) {
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
		private SortedMap<Integer, List<Polynomial<T>>> polynomialIdealGeneratorWeights;
		private PolynomialIdeal<T> asPolynomialIdeal;

		private CoordinateIdeal(CoordinateRing<T> ring, List<CoordinateRingElement<T>> generators) {
			super(ring);
			List<Polynomial<T>> generatorPolynomials = new ArrayList<>();
			for (CoordinateRingElement<T> generator : generators) {
				generatorPolynomials.add(generator.polynomial);
			}
			generatorPolynomials.addAll(ring.getIdeal().generators());
			this.ring = ring;
			PolynomialRing<T> p = ring.getPolynomialRing();
			this.asPolynomialIdeal = p.getIdeal(generatorPolynomials);
			this.generators = new ArrayList<>();
			this.polynomialIdealGeneratorWeights = new TreeMap<>();
			for (Polynomial<T> idealGenerator : ring.ideal.generators()) {
				List<Polynomial<T>> generate = asPolynomialIdeal.generate(idealGenerator);
				int index = 0;
				int last = polynomialIdealGeneratorWeights.isEmpty() ? 0
						: polynomialIdealGeneratorWeights.lastKey() + 1;
				boolean found = false;
				for (int i = last; i < generate.size(); i++) {
					if (ring.getIdeal().contains(asPolynomialIdeal.generators().get(i))) {
						continue;
					}
					index = i;
					if (p.isUnit(generate.get(index))) {
						found = true;
						break;
					}
				}
				if (found) {
					Polynomial<T> unit = p.negative(p.inverse(generate.get(index)));
					List<Polynomial<T>> coefficients = new ArrayList<>();
					for (int i = 0; i < generate.size(); i++) {
						if (index == i) {
							coefficients.add(p.zero());
							continue;
						}
						coefficients.add(p.multiply(unit, generate.get(i)));
					}
					polynomialIdealGeneratorWeights.put(index, coefficients);
				}
			}
			for (int i = 0; i < asPolynomialIdeal.generators().size(); i++) {
				if (polynomialIdealGeneratorWeights.containsKey(i)) {
					continue;
				}
				Polynomial<T> generator = asPolynomialIdeal.generators().get(i);
				Polynomial<T> reduced = ring.getIdeal().residue(generator);
				if (!reduced.equals(ring.getPolynomialRing().zero())) {
					this.generators.add(ring.getEmbedding(reduced));
				} else {
					polynomialIdealGeneratorWeights.put(i, Collections.emptyList());
				}
			}
		}

		@Override
		public boolean isPrime() {
			return asPolynomialIdeal.isPrime();
		}

		@Override
		public boolean isMaximal() {
			return asPolynomialIdeal.isMaximal();
		}

		public boolean isMaximalOverAlgebraicClosure() {
			return asPolynomialIdeal.isMaximalOverAlgebraicClosure();
		}

		public PolynomialIdeal<T> asPolynomialIdeal() {
			return asPolynomialIdeal;
		}

		@Override
		public List<CoordinateRingElement<T>> generators() {
			return Collections.unmodifiableList(generators);
		}

		@Override
		public List<CoordinateRingElement<T>> generate(CoordinateRingElement<T> t) {
			List<Polynomial<T>> generate = this.asPolynomialIdeal.generate(t.polynomial);
			PolynomialRing<T> p = ring.getPolynomialRing();
			List<Polynomial<T>> preResult = new ArrayList<>();
			for (int i = 0; i < generate.size(); i++) {
				preResult.add(generate.get(i));
			}
			boolean done = false;
			while (!done) {
				done = true;
				for (int i = 0; i < generate.size(); i++) {
					if (polynomialIdealGeneratorWeights.containsKey(i)) {
						List<Polynomial<T>> cofactors = polynomialIdealGeneratorWeights.get(i);
						Polynomial<T> value = preResult.get(i);
						preResult.set(i, p.zero());
						if (cofactors.isEmpty()) {
							continue;
						}
						if (value.equals(p.zero())) {
							continue;
						}
						for (int j = 0; j < generate.size(); j++) {
							if (i == j) {
								continue;
							}
							if (j < i && !cofactors.get(j).equals(p.zero())
									&& polynomialIdealGeneratorWeights.containsKey(j)) {
								done = false;
							}
							preResult.set(j, p.add(preResult.get(j), p.multiply(cofactors.get(j), value)));
						}
					}
				}
			}
			List<CoordinateRingElement<T>> result = new ArrayList<>();
			for (int i = 0; i < preResult.size(); i++) {
				if (!polynomialIdealGeneratorWeights.containsKey(i)) {
					result.add(ring.getEmbedding(preResult.get(i)));
				} else if (!preResult.get(i).equals(p.zero())) {
					throw new ArithmeticException("Did not clear all pre results!");
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

		public PolynomialIdeal<T> getPolynomialIdeal() {
			return asPolynomialIdeal;
		}

	}
}
