package fields.helper;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.AlgebraicRingExtension;
import fields.interfaces.AlgebraicRingExtension.PolynomialRingAsCoordinateRing;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.Vector;

public class AlgebraicRingExtensionUnivariatePolynomialRing<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends AlgebraicRingExtension<T, S, Ext>>
		extends AbstractPolynomialRing<S> implements UnivariatePolynomialRing<S> {
	private PolynomialRingAsCoordinateRing<T, S> asCoordinateRing;
	private CoordinateRing<T> coordinateRing;
	private PolynomialRing<T> multivariatePolynomialRing;
	private UnivariatePolynomialRing<T> basePolynomialRing;
	private Ext algebraicRingExtension;
	private Ring<T> base;
	private ExtensionPolynomial zero;
	private ExtensionPolynomial one;
	private ExtensionPolynomial var;
	private SortedSet<Monomial> monomials;
	

	private class ExtensionPolynomial extends AbstractElement<Polynomial<S>> implements UnivariatePolynomial<S> {
		private CoordinateRingElement<T> value;
		
		private ExtensionPolynomial(CoordinateRingElement<T> value) {
			this.value = value;
		}

		@Override
		public int compareTo(Polynomial<S> o) {
			if (!(o instanceof AlgebraicRingExtensionUnivariatePolynomialRing.ExtensionPolynomial)) {
				throw new RuntimeException("mixing polynomial types");
			}
			return value.compareTo(((ExtensionPolynomial) o).value);
		}

		@Override
		public PolynomialRing<S> getPolynomialRing() {
			return AlgebraicRingExtensionUnivariatePolynomialRing.this;
		}

		@Override
		public Ring<S> getRing() {
			return algebraicRingExtension;
		}

		@Override
		public int numberOfVariables() {
			return 1;
		}

		@Override
		public int degree() {
			return value.getElement().leadingMonomial().exponents()[0];
		}

		@Override
		public Value order() {
			for (Monomial m : value.getElement().monomials()) {
				if (!value.getElement().coefficient(m).equals(base.zero())) {
					return new Value(m.exponents()[0]);
				}
			}
			return Value.INFINITY;
		}
		
		@Override
		public int degree(int variable) {
			return degree();
		}

		@Override
		public Value order(int variable) {
			return order();
		}
		
		@Override
		public S univariateCoefficient(int degree) {
			return algebraicRingExtension.fromPolynomial(basePolynomialRing.toUnivariate(
					multivariatePolynomialRing.asUnivariatePolynomial(value.getElement(), 1).univariateCoefficient(degree)));
		}

		@Override
		public S coefficient(Monomial m) {
			return univariateCoefficient(m.exponents()[0]);
		}

		@Override
		public S leadingCoefficient() {
			return univariateCoefficient(degree());
		}

		@Override
		public Monomial leadingMonomial() {
			return getMonomial(new int[] { degree() });
		}

		@Override
		public SortedSet<Monomial> monomials() {
			return AlgebraicRingExtensionUnivariatePolynomialRing.this.monomials(degree());
		}

		@Override
		public String toString(String variable, boolean ascending) {
			if (this.value.equals(coordinateRing.zero()))
				return "0";
			StringBuffer buf = new StringBuffer();
			boolean first = true;
			for (int index = 0; index <= degree(); index++) {
				int i = index;
				if (!ascending) {
					i = degree() - index;
				}
				if (this.univariateCoefficient(i).equals(algebraicRingExtension.zero())) {
					continue;
				}
				if (first) {
					first = false;
				} else {
					buf.append(" + ");
				}
				String coefficient = this.univariateCoefficient(i).toString();
				if (coefficient.contains(" ")) {
					coefficient = "(" + coefficient + ")";
				}
				if (i == 0) {
					buf.append(coefficient);
				} else {
					if (this.univariateCoefficient(i).equals(algebraicRingExtension.one())) {
						buf.append(variable);
					} else {
						buf.append(coefficient + "*" + variable);
					}
					if (i > 1) {
						buf.append("^" + i);
					}
				}
			}
			return buf.toString();
		}

		@Override
		public String toString(String[] variables, boolean ascending) {
			return toString(variables[0], ascending);
		}

	}

	public AlgebraicRingExtensionUnivariatePolynomialRing(Ext algebraicRingExtension) {
		super(false);
		this.algebraicRingExtension = algebraicRingExtension;
		this.asCoordinateRing=algebraicRingExtension.asCoordinateRing(this);
		this.coordinateRing = asCoordinateRing.getCoordinateRing();
		this.multivariatePolynomialRing = coordinateRing.getPolynomialRing();
		this.basePolynomialRing = algebraicRingExtension.getRing().getUnivariatePolynomialRing();
		this.zero = new ExtensionPolynomial(coordinateRing.zero());
		this.one = new ExtensionPolynomial(coordinateRing.one());
		this.var = new ExtensionPolynomial(coordinateRing.getVar(1));
		this.monomials = new TreeSet<>();
	}
	
	@Override
	public Exactness exactness() {
		return algebraicRingExtension.exactness();
	}

	@Override
	public Ext getRing() {
		return algebraicRingExtension;
	}

	@Override
	public int numberOfVariables() {
		return 1;
	}

	@Override
	public Comparator<Monomial> getComparator() {
		return Monomial.LEX;
	}

	@Override
	public UnivariatePolynomial<S> zero() {
		return zero;
	}

	@Override
	public UnivariatePolynomial<S> one() {
		return one;
	}

	@Override
	public UnivariatePolynomial<S> getVar() {
		return var;
	}

	@Override
	public Polynomial<S> getVar(int i) {
		return getVar();
	}

	@Override
	public UnivariatePolynomial<S> getVarPower(int power) {
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getVarPower(1, power)));
	}

	@Override
	public UnivariatePolynomial<S> getVarPower(int i, int power) {
		return getVarPower(power);
	}

	@Override
	public UnivariatePolynomial<S> getEmbedding(S t) {
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getEmbedding(t.asPolynomial(), new int[] { 1 })));
	}

	@Override
	public UnivariatePolynomial<S> getPolynomial(List<S> coefficients) {
		SortedMap<Monomial, T> c = new TreeMap<>();
		for (int i = 0; i < coefficients.size(); i++) {
			S s = coefficients.get(i);
			UnivariatePolynomial<T> p = s.asPolynomial();
			for (int j = 0; j <= p.degree(); j++) {
				c.put(multivariatePolynomialRing.getMonomial(new int[] { i, j }), p.univariateCoefficient(j));
			}
		}
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getPolynomial(c)));
	}

	@Override
	public UnivariatePolynomial<S> getPolynomial(@SuppressWarnings("unchecked") S... coefficients) {
		return getPolynomial(Arrays.asList(coefficients));
	}

	@Override
	public UnivariatePolynomial<S> getPolynomial(Map<Monomial, S> t) {
		SortedMap<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.keySet()) {
			S s = t.get(m);
			for (int i = 0; i <= s.asPolynomial().degree(); i++) {
				c.put(multivariatePolynomialRing.getMonomial(new int[] { m.exponents()[0], i }),
						s.asPolynomial().univariateCoefficient(i));
			}
		}
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getPolynomial(c)));
	}

	@Override
	public UnivariatePolynomial<S> getEmbedding(S t, int exponent) {
		SortedMap<Monomial, T> c = new TreeMap<>();
		for (int i = 0; i <= t.asPolynomial().degree(); i++) {
			c.put(multivariatePolynomialRing.getMonomial(new int[] { exponent, i }),
					t.asPolynomial().univariateCoefficient(i));
		}
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getPolynomial(c)));
	}

	@Override
	public UnivariatePolynomial<S> getEmbedding(S t, int[] exponents) {
		return getEmbedding(t, exponents[0]);
	}

	@Override
	public UnivariatePolynomial<S> getLinear(List<S> coeff) {
		if (coeff.size() != 1)
			throw new ArithmeticException("size mismatch!");
		return this.getEmbedding(coeff.get(0), 1);
	}

	@Override
	public UnivariatePolynomial<S> toUnivariate(Polynomial<S> t) {
		if (t instanceof AlgebraicRingExtensionUnivariatePolynomialRing.ExtensionPolynomial) {
			return (ExtensionPolynomial) t;
		}
		return getEmbedding(t, new int[] { 0 });
	}

	@Override
	public boolean isHomogeneous(Polynomial<S> t) {
		if (t.equals(zero)) {
			return true;
		}
		UnivariatePolynomial<S> p = toUnivariate(t);
		return p.degree() == p.order().value();
	}

	@Override
	public PolynomialRing<S> addVariableWithElimination(int shift) {
		return AbstractPolynomialRing.getPolynomialRing(algebraicRingExtension, shift + 1,
				new Monomial.EliminationOrder(Monomial.LEX, Monomial.LEX, shift));
	}

	@Override
	public UnivariatePolynomial<S> getEmbedding(Polynomial<S> t, int[] map) {
		Map<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			int newexp = -1;
			for (int i = 0; i < t.numberOfVariables(); i++) {
				if (m.exponents()[i] != 0 && map[i] == 0) {
					newexp = m.exponents()[i];
				}
			}
			S s = t.coefficient(m);
			for (int i = 0; i <= s.asPolynomial().degree(); i++) {
				c.put(multivariatePolynomialRing.getMonomial(new int[] { newexp, i }),
						s.asPolynomial().univariateCoefficient(i));
			}
		}
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getPolynomial(c)));
	}

	@Override
	public <U extends Element<U>> UnivariatePolynomial<S> getEmbedding(Polynomial<U> t, MathMap<U, S> mapping) {
		SortedMap<Monomial, T> c = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			S s = mapping.evaluate(t.coefficient(m));
			for (int i = 0; i <= s.asPolynomial().degree(); i++) {
				c.put(getMonomial(new int[] { m.exponents()[0], i }), s.asPolynomial().univariateCoefficient(i));
			}
		}
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getPolynomial(c)));
	}

	@Override
	public Polynomial<S> getRandomElement() {
		return getRandomElement(new Random().nextInt(5));
	}

	@Override
	public UnivariatePolynomial<S> getRandomElement(int degree) {
		return new ExtensionPolynomial(coordinateRing.getEmbedding(multivariatePolynomialRing.getRandomElement(degree)));
	}

	@Override
	public SortedSet<Monomial> monomials(int degree) {
		for (int i = monomials.last().exponents()[0]; i < degree; i++) {
			monomials.add(getMonomial(new int[] { i + 1 }));
		}
		return monomials.headSet(getMonomial(new int[] { degree + 1 }));
	}

	@Override
	public Monomial getMonomial(int[] exponents) {
		return new Monomial(Monomial.LEX, exponents);
	}

	private ExtensionPolynomial toExtensionPolynomial(Polynomial<S> s) {
		if (s instanceof AlgebraicRingExtensionUnivariatePolynomialRing.ExtensionPolynomial) {
			return (ExtensionPolynomial)s;
		}
		return toExtensionPolynomial(getEmbedding(s));
	}
	
	@Override
	public UnivariatePolynomial<S> add(Polynomial<S> s1, Polynomial<S> s2) {
		return new ExtensionPolynomial(coordinateRing.add(toExtensionPolynomial(s1).value, toExtensionPolynomial(s2).value));
	}

	@Override
	public UnivariatePolynomial<S> negative(Polynomial<S> s) {
		return new ExtensionPolynomial(coordinateRing.negative(toExtensionPolynomial(s).value));
	}

	@Override
	public UnivariatePolynomial<S> multiply(Polynomial<S> t1, Polynomial<S> t2) {
		return new ExtensionPolynomial(coordinateRing.multiply(toExtensionPolynomial(t1).value, toExtensionPolynomial(t2).value));
	}

	@Override
	public UnivariatePolynomial<S> multiplyPower(int power, Polynomial<S> t) {
		return multiply(getVarPower(power), t);
	}

	@Override
	public UnivariatePolynomial<S> multiply(S t1, Polynomial<S> t2) {
		return multiply(getEmbedding(t1), t2);
	}

	@Override
	public UnivariatePolynomial<S> multiplyPower(int power, S t1, Polynomial<S> t2) {
		return multiply(multiplyPower(power, getEmbedding(t1)), t2);
	}

	@Override
	public UnivariatePolynomial<S> divideScalar(Polynomial<S> dividend, S divisor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isZeroDivisor(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isUnit(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public UnivariatePolynomial<S> inverse(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isDivisible(Polynomial<S> dividend, Polynomial<S> divisor) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public GeneralQuotientAndRemainderResult<S> generalQuotientAndRemainder(Polynomial<S> dividend,
			List<Polynomial<S>> divisors) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> pseudoRemainder(Polynomial<S> dividend, Polynomial<S> divisor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigInteger euclidMeasure(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<Polynomial<S>> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public FactorizationResult<Polynomial<S>, S> squareFreeFactorization(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public List<Polynomial<S>> getAlgebraGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<S> asVector(Polynomial<S> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<S> asVector(Polynomial<S> s, int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> fromVector(Vector<S> vector) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<Polynomial<S>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<List<S>> nonTrivialCombinations(List<Polynomial<S>> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<Polynomial<S>> iterator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<UnivariatePolynomial<S>> polynomials(int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<UnivariatePolynomial<S>> polynomialSet(int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<UnivariatePolynomial<S>> monicPolynomials(int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<UnivariatePolynomial<S>> monicPolynomialSet(int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public S evaluate(Polynomial<S> t, List<S> ts) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Polynomial<S> partiallyEvaluate(Polynomial<S> t, List<S> ts) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> substitute(Polynomial<S> t, List<Polynomial<S>> values) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean canBeNormalized(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public UnivariatePolynomial<S> normalize(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public S content(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> contentFree(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> derivative(Polynomial<S> t, int variable) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> derivative(Polynomial<S> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<S> trace(Polynomial<S> t, int degree) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public S resultant(Polynomial<S> t1, Polynomial<S> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Polynomial<S> resultant(Polynomial<S> t1, Polynomial<S> t2, int variable) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public PolynomialRing<S> eliminateVariable() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<Polynomial<S>> asUnivariatePolynomial(Polynomial<S> t, int variable) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Polynomial<S> fromUnivariatePolynomial(UnivariatePolynomial<Polynomial<S>> t, int variable) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> gcd(Polynomial<S> t1, Polynomial<S> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public S discriminant(Polynomial<S> t1) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ExtendedResultantResult<S> extendedResultant(Polynomial<S> t1, Polynomial<S> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public UnivariatePolynomial<S> round(Polynomial<S> t, int degree) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Vector<S> asVector(Polynomial<S> t, int[] degrees) {
		// TODO Auto-generated method stub
		return null;
	}
	
	@Override
	public Polynomial<S> fromVector(Vector<S> t, int[] degrees) {
		// TODO Auto-generated method stub
		return null;
	}
}
