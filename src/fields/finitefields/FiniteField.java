package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PrimeFieldElement;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.helper.ExtensionField;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.UnivariatePolynomialRing;
import util.Pair;

public class FiniteField extends AbstractField<FFE> {
	private PrimeField base;
	private ExtensionField<PrimeFieldElement> field;
	private Map<Polynomial<FFE>, List<Pair<Polynomial<FFE>, Integer>>> distinctDegreeFactorizatonCache;
	private Map<Polynomial<FFE>, Integer> distinctDegreeFactorizationMaxDegree;
	private Map<Polynomial<FFE>, Map<Integer, List<Polynomial<FFE>>>> irreducibleFactorizatonCache;

	public static class FFE extends AbstractElement<FFE> {
		private ExtensionFieldElement<PrimeFieldElement> e;

		private FFE(ExtensionFieldElement<PrimeFieldElement> e) {
			this.e = e;
		}

		@Override
		public int compareTo(FFE o) {
			return e.compareTo(o.e);
		}
		
		public ExtensionFieldElement<PrimeFieldElement> asExtensionFieldElement() {
			return e;
		}
		
		public static FFE fromExtensionFieldElement(ExtensionFieldElement<PrimeFieldElement> t) {
			return new FFE(t);
		}

		@Override
		public String toString() {
			return e.toString();
		}
	}

	private static Polynomial<PrimeFieldElement> findIrreduciblePolynomial(int degree, PrimeField field) {
		Polynomial<PrimeFieldElement> result;
		PolynomialRing<PrimeFieldElement> r = field.getUnivariatePolynomialRing();
		do {
			Map<Monomial, PrimeFieldElement> coeffs = new TreeMap<>();
			coeffs.put(r.getMonomial(new int[] { degree }), field.one());
			for (int i = 0; i < degree; i++) {
				coeffs.put(r.getMonomial(new int[] { i }), field.getRandomElement());
			}
			result = r.getPolynomial(coeffs);
		} while (field.factorization(result).size() != 1);
		return result;
	}

	public FiniteField(PrimeField base) {
		this(base.getUnivariatePolynomialRing().getVar(1), base);
	}

	public FiniteField(PrimeField base, int degree) {
		this(findIrreduciblePolynomial(degree, base), base);
	}

	public FiniteField(Polynomial<PrimeFieldElement> minimalpolynomial, PrimeField base) {
		this.base = base;
		this.field = new ExtensionField<>(minimalpolynomial, base);
		this.distinctDegreeFactorizatonCache = new TreeMap<>();
		this.distinctDegreeFactorizationMaxDegree = new TreeMap<>();
		this.irreducibleFactorizatonCache = new TreeMap<>();
	}
	
	public Polynomial<PrimeFieldElement> minimalPolynomial() {
		return field.minimalPolynomial();
	}
	
	public ExtensionField<PrimeFieldElement> asExtensionField() {
		return field;
	}

	@Override
	public FFE zero() {
		return new FFE(field.zero());
	}

	@Override
	public FFE one() {
		return new FFE(field.one());
	}

	public FFE getEmbedding(PrimeFieldElement t) {
		return new FFE(field.getEmbedding(t));
	}

	@Override
	public BigInteger characteristic() {
		return base.characteristic();
	}

	@Override
	public FFE add(FFE t1, FFE t2) {
		return new FFE(field.add(t1.e, t2.e));
	}

	@Override
	public FFE negative(FFE t) {
		return new FFE(field.negative(t.e));
	}

	@Override
	public FFE multiply(FFE t1, FFE t2) {
		return new FFE(field.multiply(t1.e, t2.e));
	}

	@Override
	public FFE inverse(FFE t) {
		return new FFE(field.inverse(t.e));
	}

	@Override
	public FFE getRandomElement() {
		return new FFE(field.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return field.getNumberOfElements();
	}

	@Override
	public Iterator<FFE> iterator() {
		return new Iterator<FFE>() {
			Iterator<ExtensionFieldElement<PrimeFieldElement>> it = field.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public FFE next() {
				return new FFE(it.next());
			}
		};
	}

	public List<Polynomial<FFE>> factorization(Polynomial<FFE> t) {
		UnivariatePolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Multivariate factorization not implemented.");
		}
		List<Polynomial<FFE>> result = new ArrayList<>();
		if (!t.leadingCoefficient().equals(this.one())) {
			result.add(ring.getEmbedding(t.leadingCoefficient()));
			t = ring.normalize(t);
		}
		List<Polynomial<FFE>> squareFreeFactors = ring.squareFreeFactorization(t);
		for (Polynomial<FFE> sff : squareFreeFactors) {
			result.addAll(factorizeSquareFree(sff));
		}
		return result;
	}

	private List<Polynomial<FFE>> factorizeSquareFree(Polynomial<FFE> t) {
		List<Polynomial<FFE>> result = new ArrayList<>();
		List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(t, -1);
		for (Pair<Polynomial<FFE>, Integer> ddf : distinctDegreeFactors) {
			result.addAll(this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond()));
		}
		return result;
	}

	/**
	 * Factorizes a square free polynomial into factors that of a distinct degree.
	 * 
	 * @param t         univariate, normalized, square free polynomial over a finite
	 *                  field.
	 * @param maxDegree maximum degree of results needed. Set to -1 for all factors.
	 * @return list of pairs of polynomials and integers, such that the polynomials
	 *         are factors of t, that each consists of irreducible factors of the
	 *         the given degree.
	 */
	public List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactorization(Polynomial<FFE> t, int maxDegree) {
		if (this.distinctDegreeFactorizatonCache.containsKey(t)) {
			int previousMaxDegree = this.distinctDegreeFactorizationMaxDegree.get(t);
			if (previousMaxDegree == -1 || previousMaxDegree >= maxDegree) {
				return this.distinctDegreeFactorizatonCache.get(t);
			}
		}
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!t.leadingCoefficient().equals(this.one())) {
			throw new RuntimeException("Not normalized.");
		}
		PolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		List<Pair<Polynomial<FFE>, Integer>> result = new ArrayList<>();
		Polynomial<FFE> work = t;
		CoordinateRing<FFE> cr = new CoordinateRing<>(ring, ring.getIdeal(Collections.singletonList(work)));
		CoordinateRingElement<FFE> x = cr.getEmbedding(ring.getVar(1));
		CoordinateRingElement<FFE> xqi = x;
		for (int degree = 1; work.degree() >= 2 * degree; degree++) {
			if (maxDegree >= 0 && degree > maxDegree) {
				break;
			}
			xqi = cr.power(xqi, this.getNumberOfElements());
			CoordinateRingElement<FFE> xqimx = cr.subtract(xqi, x);
			Polynomial<FFE> gcd = ring.normalize(ring.gcd(work, xqimx.getElement()));
			if (!gcd.equals(ring.one())) {
				result.add(new Pair<>(gcd, degree));
				work = ring.divide(work, gcd);
				cr = new CoordinateRing<FFE>(ring, ring.getIdeal(Collections.singletonList(work)));
				x = cr.getEmbedding(ring.getVar(1));
				xqi = x;
				for (int j = 0; j < degree; j++) {
					xqi = cr.power(xqi, this.getNumberOfElements());
				}
			}
		}
		if (!work.equals(ring.one())) {
			result.add(new Pair<>(work, work.degree()));
		}
		this.distinctDegreeFactorizatonCache.put(t, result);
		this.distinctDegreeFactorizationMaxDegree.put(t, maxDegree);
		return result;
	}

	private List<Polynomial<FFE>> irreducibleFactorization(Polynomial<FFE> t, int degree) {
		if (!this.irreducibleFactorizatonCache.containsKey(t)) {
			this.irreducibleFactorizatonCache.put(t, new TreeMap<>());
		}
		Map<Integer, List<Polynomial<FFE>>> cache = this.irreducibleFactorizatonCache.get(t);
		if (cache.containsKey(degree)) {
			return cache.get(degree);
		}
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not univariate.");
		}
		if (!t.leadingCoefficient().equals(this.one())) {
			throw new RuntimeException("Not normalized.");
		}
		if (t.degree() % degree != 0) {
			throw new RuntimeException("Inputs impossible!");
		}
		PolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		int numberOfFactors = t.degree() / degree;
		CoordinateRing<FFE> cr = new CoordinateRing<>(ring, ring.getIdeal(Collections.singletonList(t)));
		List<Polynomial<FFE>> factors = new ArrayList<>();
		factors.add(t);
		while (factors.size() < numberOfFactors) {
			// Generate random polynomial with degree smaller than t.
			Polynomial<FFE> h = ring.zero();
			for (int i = 0; i < t.degree(); i++) {
				h = ring.add(h, ring.multiply(this.getRandomElement(), ring.getVarPower(1, i)));
			}
			BigInteger one = BigInteger.ONE;
			BigInteger exponent = this.getNumberOfElements().pow(degree).subtract(one).shiftRight(1);
			Polynomial<FFE> g = cr.subtract(cr.power(cr.getEmbedding(h), exponent), cr.one()).getElement();
			List<Polynomial<FFE>> newFactors = new ArrayList<Polynomial<FFE>>();
			for (Polynomial<FFE> factor : factors) {
				if (factor.degree() == degree) {
					newFactors.add(factor);
					continue;
				}
				Polynomial<FFE> gcd = ring.normalize(ring.gcd(factor, g));
				if (!gcd.equals(ring.one()) && !gcd.equals(factor)) {
					newFactors.add(gcd);
					newFactors.add(ring.quotientAndRemainder(factor, gcd).get(0));
				} else {
					newFactors.add(factor);
				}
			}
			factors = newFactors;
		}
		cache.put(degree, factors);
		return factors;
	}

	@Override
	public boolean hasRoots(Polynomial<FFE> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Not implemented.");
		}
		PolynomialRing<FFE> ring = this.getUnivariatePolynomialRing();
		CoordinateRing<FFE> cr = new CoordinateRing<FFE>(ring, ring.getIdeal(Collections.singletonList(t)));
		CoordinateRingElement<FFE> x = cr.getEmbedding(ring.getVar(1));
		CoordinateRingElement<FFE> xq = cr.power(x, this.getNumberOfElements());
		Polynomial<FFE> xqmx = cr.subtract(xq, x).getElement();
		List<Polynomial<FFE>> egcd = ring.extendedEuclidean(t, xqmx);
		return egcd.get(0).degree() > 0;
	}
	
	@Override
	public List<FFE> roots(Polynomial<FFE> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Multivariate root finding not implemented. Use ideal logic.");
		}
		UnivariatePolynomialRing<FFE> ring = getUnivariatePolynomialRing();
		List<FFE> result = new ArrayList<FFE>();
		t = ring.normalize(t);
		List<Polynomial<FFE>> squareFreeFactors = ring.squareFreeFactorization(t);
		for (Polynomial<FFE> sff : squareFreeFactors) {
			List<Pair<Polynomial<FFE>, Integer>> distinctDegreeFactors = this.distinctDegreeFactorization(sff, 1);
			for (Pair<Polynomial<FFE>, Integer> ddf : distinctDegreeFactors) {
				if (ddf.getSecond() != 1) {
					continue;
				}
				List<Polynomial<FFE>> linear = this.irreducibleFactorization(ddf.getFirst(), ddf.getSecond());
				for (Polynomial<FFE> linearFactor : linear) {
					result.add(this.negative(linearFactor.coefficient(ring.getMonomial(new int[] { 0 }))));
				}
			}
		}
		return result;
	}
	
	@Override
	public String toString() {
		return field.toString();
	}
}
