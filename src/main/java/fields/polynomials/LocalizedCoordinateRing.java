package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.FieldOfFractions;
import fields.helper.FieldOfFractions.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.local.FormalPowerSeries;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.local.LocalRingImplementation;
import fields.local.Value;
import fields.local.ValueGroup;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import util.Identity;
import util.Pair;

public class LocalizedCoordinateRing<T extends Element<T>> extends AbstractField<LocalizedElement<T>>
		implements DiscreteValuationField<LocalizedElement<T>, T> {
	public static class LocalizedElement<T extends Element<T>> extends AbstractElement<LocalizedElement<T>> {
		private Fraction<Polynomial<T>> asPolynomialFraction;
		private boolean canonical;
		private Value value;
		private Fraction<Polynomial<T>> normalized;
		private CoordinateRing<T> coordinateRing;
		private LocalizedCoordinateRing<T> localizedRing;

		private LocalizedElement(LocalizedCoordinateRing<T> localizedRing, Fraction<Polynomial<T>> asPolynomialFraction,
				CoordinateRing<T> coordinateRing, FieldOfFractions<Polynomial<T>> fieldOfFractions) {
			this.coordinateRing = coordinateRing;
			this.asPolynomialFraction = asPolynomialFraction;
			this.localizedRing = localizedRing;
			this.canonical = false;
		}

		public void canonicalize() {
			if (canonical) {
				return;
			}
			Pair<Fraction<Polynomial<T>>, Value> normalized = localizedRing.asUniformizerPower(this);
			this.value = normalized.getSecond();
			this.normalized = normalized.getFirst();
			canonical = true;
		}

		public CoordinateRingElement<T> getNumerator() {
			return coordinateRing.getEmbedding(asPolynomialFraction.getNumerator());
		}

		public CoordinateRingElement<T> getDenominator() {
			return coordinateRing.getEmbedding(asPolynomialFraction.getDenominator());
		}

		@Override
		public boolean equals(Object o) {
			if (!(o instanceof LocalizedElement<?>)) {
				return false;
			}
			@SuppressWarnings("unchecked")
			LocalizedElement<T> other = (LocalizedElement<T>) o;
			PolynomialRing<T> p = coordinateRing.getPolynomialRing();
			Polynomial<T> lhs = p.multiply(asPolynomialFraction.getNumerator(),
					other.asPolynomialFraction.getDenominator());
			Polynomial<T> rhs = p.multiply(other.asPolynomialFraction.getNumerator(),
					asPolynomialFraction.getDenominator());
			return coordinateRing.getEmbedding(lhs).equals(coordinateRing.getEmbedding(rhs));
		}

		@Override
		public int compareTo(LocalizedElement<T> o) {
			canonicalize();
			o.canonicalize();
			int cmp = value.compareTo(o.value);
			if (cmp != 0) {
				return cmp;
			}
			asPolynomialFraction.canonicalize();
			o.asPolynomialFraction.canonicalize();
			PolynomialRing<T> p = coordinateRing.getPolynomialRing();
			Polynomial<T> lhs = p.multiply(normalized.getNumerator(), o.normalized.getDenominator());
			Polynomial<T> rhs = p.multiply(o.normalized.getNumerator(), normalized.getDenominator());
			return coordinateRing.getEmbedding(lhs).compareTo(coordinateRing.getEmbedding(rhs));
		}

		@Override
		public String toString() {
			return asPolynomialFraction.toString();
		}

		public Fraction<Polynomial<T>> asPolynomialFraction() {
			return asPolynomialFraction;
		}

	}

	private FieldOfFractions<Polynomial<T>> fieldOfFractions;
	private Field<T> field;
	private CoordinateRing<T> ring;
	private CoordinateIdeal<T> ideal;
	private DiscreteValuationRing<LocalizedElement<T>, T> localRing;
	private LocalizedElement<T> uniformizer;
	private Polynomial<T> uniformizerAsPolynomial;
	private int uniformizerVariable;
	private PolynomialRing<T> polynomialRing;
	private Reals r = Reals.r(1024);

	public LocalizedCoordinateRing(Field<T> field, CoordinateRing<T> ring, CoordinateIdeal<T> ideal) {
		this.field = field;
		PolynomialRing<T> polynomialRing = ring.getPolynomialRing();
		CoordinateRing<T> reduction = ideal.divideOut();
		if (ring.krullDimension() != 1) {
			throw new ArithmeticException("not a curve");
		}
		if (reduction.krullDimension() != ring.krullDimension() - 1) {
			throw new ArithmeticException("wrong dimensions");
		}
		if (!ideal.isMaximalOverAlgebraicClosure()) {
			throw new ArithmeticException("not a maximal ideal!");
		}
		this.uniformizerVariable = -1;
		List<List<T>> differentialFormEquations = new ArrayList<>();
		for (Polynomial<T> generator : ring.getIdeal().generators()) {
			List<T> differentialFormEquation = new ArrayList<>();
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				Polynomial<T> residue = ideal.asPolynomialIdeal().residue(polynomialRing.derivative(generator, i + 1));
				if (residue.degree() > 0) {
					throw new ArithmeticException("ideal not maximal!");
				}
				differentialFormEquation.add(residue.leadingCoefficient());
			}
			differentialFormEquations.add(differentialFormEquation);
		}
		if (!differentialFormEquations.isEmpty()) {
			Matrix<T> differentialForms = new Matrix<>(differentialFormEquations);
			MatrixModule<T> mm = new MatrixModule<>(field, ring.getIdeal().generators().size(),
					polynomialRing.numberOfVariables());
			FiniteVectorSpace<T> space = new FiniteVectorSpace<>(field, polynomialRing.numberOfVariables());
			FiniteVectorSpace<T> space2 = new FiniteVectorSpace<>(field, ring.getIdeal().generators().size());
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				if (mm.multiply(differentialForms, space.getUnitVector(i + 1)).equals(space2.zero())) {
					this.uniformizerVariable = i + 1;
					break;
				}
			}
		}
		if (this.uniformizerVariable == -1) {
			this.uniformizerVariable = 1;
		}
		this.polynomialRing = AbstractPolynomialRing.getPolynomialRing(field, polynomialRing.numberOfVariables(),
				new Monomial.InvertedEliminateVariableOrder(uniformizerVariable));
		this.fieldOfFractions = new FieldOfFractions<>(this.polynomialRing);
		List<Polynomial<T>> idealGenerators = new ArrayList<>();
		for (Polynomial<T> generator : ring.getIdeal().generators()) {
			idealGenerators.add(this.polynomialRing.getEmbedding(generator));
		}
		this.ring = this.polynomialRing.getIdeal(idealGenerators).divideOut();
		List<CoordinateRingElement<T>> generators = new ArrayList<>();
		for (Polynomial<T> generator : ideal.asPolynomialIdeal().generators()) {
			generators.add(this.ring.getEmbedding(this.polynomialRing.getEmbedding(generator)));
		}
		this.ideal = this.ring.getIdeal(generators);
		this.uniformizerAsPolynomial = this.polynomialRing.subtract(this.polynomialRing.getVar(uniformizerVariable),
				this.ideal.asPolynomialIdeal().residue(this.polynomialRing.getVar(uniformizerVariable)));
		this.uniformizer = getEmbedding(uniformizerAsPolynomial);
		this.localRing = new LocalRingImplementation<>(this, ring + "_" + ideal);
	}

	private Pair<Fraction<Polynomial<T>>, Value> asUniformizerPower(LocalizedElement<T> t) {
		if (t.canonical) {
			return new Pair<>(t.normalized, t.value);
		}
		if (t.equals(zero())) {
			return new Pair<>(fieldOfFractions.zero(), Value.INFINITY);
		}
		ValueGroup g = ValueGroup.g();
		if (ideal.contains(t.getDenominator())) {
			Pair<Fraction<Polynomial<T>>, Value> numeratorPair = asUniformizerPower(getEmbedding(t.getNumerator()));
			Pair<Fraction<Polynomial<T>>, Value> denominatorPair = asUniformizerPower(getEmbedding(t.getDenominator()));
			return new Pair<>(fieldOfFractions.divide(numeratorPair.getFirst(), denominatorPair.getFirst()),
					g.subtract(numeratorPair.getSecond(), denominatorPair.getSecond()));
		}
		int power = 0;
		Polynomial<T> numerator = t.getNumerator().getElement();
		while (polynomialRing.isDivisible(numerator, uniformizerAsPolynomial)) {
			numerator = polynomialRing.divide(numerator, uniformizerAsPolynomial);
			power++;
		}
		if (!ideal.contains(ring.getEmbedding(numerator))) {
			return new Pair<>(fieldOfFractions.getFraction(numerator, t.asPolynomialFraction.getDenominator()),
					new Value(power));
		}
		List<Polynomial<T>> list = new ArrayList<>();
		list.add(numerator);
		list.add(uniformizerAsPolynomial);
		IdealResult<Polynomial<T>, PolynomialIdeal<T>> asIdeal = polynomialRing.getIdealWithTransforms(list);
		for (Polynomial<T> generator : ring.getIdeal().generators()) {
			Polynomial<T> gen = generator;
			while (!asIdeal.getIdeal().contains(generator)) {
				generator = polynomialRing.multiply(generator, gen);
			}
			if (asIdeal.getIdeal().contains(generator)) {
				List<Polynomial<T>> generate = asIdeal.getIdeal().generate(generator);
				List<Polynomial<T>> realGenerate = new ArrayList<>();
				realGenerate.add(polynomialRing.zero());
				realGenerate.add(polynomialRing.zero());
				for (int i = 0; i < generate.size(); i++) {
					realGenerate.set(0, polynomialRing.add(realGenerate.get(0),
							polynomialRing.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(0))));
					realGenerate.set(1, polynomialRing.add(realGenerate.get(1),
							polynomialRing.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(1))));
				}
				if (!realGenerate.get(0).equals(polynomialRing.zero())
						&& !realGenerate.get(1).equals(polynomialRing.zero())) {
					LocalizedElement<T> newElement = getElement(
							fieldOfFractions.getFraction(realGenerate.get(1), realGenerate.get(0)));
					Pair<Fraction<Polynomial<T>>, Value> newResult = asUniformizerPower(newElement);
					return new Pair<>(
							fieldOfFractions.divide(newResult.getFirst(),
									fieldOfFractions.getEmbedding(t.asPolynomialFraction.getDenominator())),
							g.add(new Value(power + 1), newResult.getSecond()));
				}
			}
		}
		return null;
	}

	private LocalizedElement<T> getElement(Fraction<Polynomial<T>> fraction) {
		return new LocalizedElement<>(this, fraction, ring, fieldOfFractions);
	}

	public LocalizedElement<T> getEmbedding(T t) {
		return getEmbedding(polynomialRing.getEmbedding(t));
	}

	public LocalizedElement<T> getEmbedding(CoordinateRingElement<T> t) {
		return getEmbedding(t.getElement());
	}

	public LocalizedElement<T> getEmbedding(Polynomial<T> t) {
		return getElement(fieldOfFractions.getEmbedding(polynomialRing.getEmbedding(t)));
	}

	public LocalizedElement<T> getEmbedding(Polynomial<T> numerator, Polynomial<T> denominator) {
		return getElement(fieldOfFractions.getFraction(numerator, denominator));
	}

	@Override
	public Real value(LocalizedElement<T> t) {
		if (t.equals(zero())) {
			return r.zero();
		}
		return r.exp(r.getInteger(-valuation(t).value()));
	}
	
	@Override
	public Reals getReals() {
		return r;
	}

	@Override
	public LocalizedElement<T> zero() {
		return getElement(fieldOfFractions.zero());
	}

	@Override
	public LocalizedElement<T> one() {
		return getElement(fieldOfFractions.one());
	}

	@Override
	public BigInteger characteristic() {
		return field.characteristic();
	}

	@Override
	public LocalizedElement<T> add(LocalizedElement<T> t1, LocalizedElement<T> t2) {
		return getElement(fieldOfFractions.add(t1.asPolynomialFraction, t2.asPolynomialFraction));
	}

	@Override
	public LocalizedElement<T> negative(LocalizedElement<T> t) {
		return getElement(fieldOfFractions.negative(t.asPolynomialFraction));
	}

	@Override
	public LocalizedElement<T> multiply(LocalizedElement<T> t1, LocalizedElement<T> t2) {
		return getElement(fieldOfFractions.multiply(t1.asPolynomialFraction, t2.asPolynomialFraction));
	}

	public LocalizedElement<T> multiply(T t1, LocalizedElement<T> t2) {
		return getEmbedding(polynomialRing.multiply(t1, t2.asPolynomialFraction.getNumerator()),
				t2.asPolynomialFraction.getDenominator());
	}

	@Override
	public LocalizedElement<T> inverse(LocalizedElement<T> t) {
		return getElement(fieldOfFractions.inverse(t.asPolynomialFraction));
	}

	@Override
	public Exactness exactness() {
		return ring.exactness();
	}

	@Override
	public LocalizedElement<T> getRandomElement() {
		return getElement(fieldOfFractions.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<LocalizedElement<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public boolean isComplete() {
		return false;
	}

	@Override
	public LocalizedElement<T> inverse(LocalizedElement<T> t, int accuracy) {
		return getElement(fieldOfFractions.inverse(t.asPolynomialFraction));
	}

	@Override
	public LocalizedElement<T> negative(LocalizedElement<T> t, int accuracy) {
		return getElement(fieldOfFractions.negative(t.asPolynomialFraction));
	}

	@Override
	public boolean isInteger(LocalizedElement<T> t) {
		return !ideal.contains(t.getDenominator());
	}

	@Override
	public Value valuation(LocalizedElement<T> t) {
		t.canonicalize();
		return t.value;
	}

	@Override
	public LocalizedElement<T> upToUniformizerPower(LocalizedElement<T> t) {
		t.canonicalize();
		return getElement(t.normalized);
	}

	@Override
	public LocalizedElement<T> uniformizer() {
		return uniformizer;
	}

	@Override
	public Field<T> residueField() {
		return field;
	}

	@Override
	public T reduceInteger(LocalizedElement<T> t) {
		t.canonicalize();
		if (t.value.compareTo(Value.ZERO) > 0) {
			return field.zero();
		}
		if (t.value.compareTo(Value.ZERO) < 0) {
			throw new ArithmeticException("Not an integer!");
		}
		CoordinateRingElement<T> numerator = ideal.residue(ring.getEmbedding(t.normalized.getNumerator()));
		CoordinateRingElement<T> denominator = ideal.residue(ring.getEmbedding(t.normalized.getDenominator()));
		if (numerator.getElement().degree() > 0 || denominator.getElement().degree() != 0) {
			throw new ArithmeticException("Not an integer (Oh no)!");
		}
		return field.divide(numerator.getElement().leadingCoefficient(), denominator.getElement().leadingCoefficient());
	}

	@Override
	public LocalizedElement<T> liftToInteger(T s) {
		return getEmbedding((ring.getPolynomialRing().getEmbedding(s)));
	}

	@Override
	public LocalizedElement<T> round(LocalizedElement<T> t, int accuracy) {
		if (t.equals(zero())) {
			return zero();
		}
		int minPower = valuation(t).value();
		LocalizedElement<T> uniformizerPower = power(uniformizer(), minPower);
		t = divide(t, uniformizerPower);
		LocalizedElement<T> result = zero();
		for (int i = minPower; i < accuracy; i++) {
			if (t.equals(zero())) {
				break;
			}
			T digit = reduceInteger(t);
			LocalizedElement<T> adjustedDigit = multiply(digit, uniformizerPower);
			result = add(result, adjustedDigit);
			uniformizerPower = multiply(uniformizerPower, uniformizer);
			t = divide(subtract(t, getEmbedding(digit)), uniformizer);
		}
		return result;
	}

	@Override
	public int getAccuracy() {
		throw new InfinityException();
	}

	@Override
	public DiscreteValuationField<LocalizedElement<T>, T> withAccuracy(int accuracy) {
		return this;
	}

	@Override
	public OtherVersion<LocalizedElement<T>, LocalizedElement<T>, T, LocalizedCoordinateRing<T>> exact() {
		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
	}

	@Override
	public OtherVersion<LocalizedElement<T>, PowerSeries<T>, T, FormalPowerSeries<T>> complete(int accuracy) {
		FormalPowerSeries<T> powerSeries = new FormalPowerSeries<>(field, accuracy);
		MathMap<LocalizedElement<T>, PowerSeries<T>> embedding = new MathMap<>() {

			@Override
			public PowerSeries<T> evaluate(LocalizedElement<T> t) {
				if (t.equals(zero())) {
					return powerSeries.zero();
				}
				List<T> coefficients = new ArrayList<>();
				int minPower = valuation(t).value();
				t = multiply(t, power(uniformizer(), -minPower));
				for (int i = minPower; i < accuracy; i++) {
					T c = reduceInteger(t);
					coefficients.add(c);
					t = divide(subtract(t,getEmbedding(c)), uniformizer());
					if (t.equals(zero())) {
						break;
					}
				}
				return powerSeries.getElement(field.getUnivariatePolynomialRing().getPolynomial(coefficients),
						minPower);
			}
		};
		MathMap<PowerSeries<T>, LocalizedElement<T>> rounded = new MathMap<>() {
			@Override
			public LocalizedElement<T> evaluate(PowerSeries<T> t) {
				Value value = powerSeries.valuation(t);
				if (value.isInfinite()) {
					return zero();
				}
				int power = Math.min(0, -value.value());
				Polynomial<T> denominator = polynomialRing.power(uniformizerAsPolynomial, power);
				t = powerSeries.multiply(t, powerSeries.power(powerSeries.uniformizer(), power));
				Polynomial<T> numerator = polynomialRing.substitute(powerSeries.roundToPolynomial(t, accuracy),
						Collections.singletonList(uniformizerAsPolynomial));
				return getEmbedding(numerator, denominator);
			}
		};
		return new OtherVersion<>(powerSeries, embedding, rounded);
	}

	@Override
	public LocalizedElement<T> getRandomInteger() {
		return getEmbedding(ring.getRandomElement());
	}

	@Override
	public DiscreteValuationRing<LocalizedElement<T>, T> ringOfIntegers() {
		return localRing;
	}
	
	@Override
	public DiscreteValuationFieldExtension<LocalizedElement<T>, T, ?, ?, ?, ?> getUniqueExtension(
			UnivariatePolynomial<LocalizedElement<T>> minimalPolynomial) {
	throw new UnsupportedOperationException("Not implemented!");
	}

}
