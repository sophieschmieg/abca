package varieties;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.FunctionMathMap;
import varieties.SimpleFunctionField.SimpleRationalFunction;
import varieties.affine.AffineScheme;
import varieties.projective.ProjectiveScheme;

public class SimpleFunctionField<T extends Element<T>>
		extends AbstractFieldExtension<TExt<T>, SimpleRationalFunction<T>, SimpleFunctionField<T>> {
	public static class SimpleRationalFunction<T extends Element<T>> extends AbstractElement<SimpleRationalFunction<T>>
			implements AlgebraicExtensionElement<TExt<T>, SimpleRationalFunction<T>> {
		private SimpleFunctionField<T> functionField;
		private UnivariatePolynomial<TExt<T>> asPolynomial;
		private Polynomial<T> lcm = null;
		private Polynomial<T> numerator = null;
		private Polynomial<T> denominator = null;

		private SimpleRationalFunction(SimpleFunctionField<T> functionField,
				UnivariatePolynomial<TExt<T>> asPolynomial) {
			this.functionField = functionField;
			this.asPolynomial = asPolynomial;
		}

		@Override
		public int compareTo(SimpleRationalFunction<T> o) {
			return asPolynomial.compareTo(o.asPolynomial);
		}

		@Override
		public UnivariatePolynomial<TExt<T>> asPolynomial() {
			return asPolynomial;
		}

		@Override
		public String toString() {
			return asPolynomial.toString("α", true);
		}

		public Polynomial<T> getDenominatorLcm() {
			if (lcm == null) {
				PolynomialRing<T> polynomialRing = functionField.getBaseField().polynomialRing();
				lcm = polynomialRing.one();
				for (int i = 0; i <= asPolynomial.degree(); i++) {
					lcm = polynomialRing.lcm(lcm, asPolynomial.univariateCoefficient(i).getDenominator());
				}
			}
			return lcm;
		}

		public Polynomial<T> getNumerator() {
			if (numerator != null) {
				Polynomial<T> denominatorLcm = getDenominatorLcm();
				PolynomialRing<T> polynomialRing = functionField.getBaseField().polynomialRing();
				numerator = polynomialRing.zero();
				for (int i = 0; i <= asPolynomial.degree(); i++) {
					numerator = polynomialRing.add(numerator,
							polynomialRing.multiply(asPolynomial.univariateCoefficient(i).getNumerator(),
									polynomialRing.divideChecked(denominatorLcm,
											asPolynomial.univariateCoefficient(i).getDenominator())));
				}
			}
			return numerator;
		}

		public Polynomial<T> getDenominator() {
			if (denominator != null) {
				denominator = getDenominatorLcm();
			}
			return denominator;
		}
	}

	private TranscendentalFieldExtension<T> transcendentalExtension;

//	public static <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>> SimpleFunctionFieldFromCoordinateRingOverExtension<B, E> forProjectiveVarietyOverExtensionField(
//			Ext extension, ProjectiveScheme<E> domain) {
//		int affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
//		AffineScheme<E> affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
//		return fromCoordinateRingOverExtensionField(extension, affineSlice.getCoordinateRing());
//	}

	public static <T extends Element<T>> SimpleFunctionFieldFromCoordinateRing<T> forProjectiveVariety(
			ProjectiveScheme<T> domain) {
		int affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
		AffineScheme<T> affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
		return fromCoordinateRing(affineSlice.getCoordinateRing());
	}
	/*
	 * public static class SimpleFunctionFieldFromCoordinateRingOverExtension<B
	 * extends Element<B>, E extends AlgebraicExtensionElement<B, E>> { private
	 * SimpleFunctionField<B> simpleFunctionField; private
	 * MathMap<RationalFunction<E>, SimpleRationalFunction<B>> isomorphism; private
	 * MathMap<SimpleRationalFunction<B>, RationalFunction<E>> inverseIsomorphism;
	 * 
	 * private
	 * SimpleFunctionFieldFromCoordinateRingOverExtension(SimpleFunctionField<B>
	 * simpleFunctionField, MathMap<RationalFunction<E>, SimpleRationalFunction<B>>
	 * isomorphism, MathMap<SimpleRationalFunction<B>, RationalFunction<E>>
	 * inverseIsomorphism) { this.simpleFunctionField = simpleFunctionField;
	 * this.isomorphism = isomorphism; this.inverseIsomorphism = inverseIsomorphism;
	 * }
	 * 
	 * public SimpleFunctionField<B> getSimpleFunctionField() { return
	 * simpleFunctionField; }
	 * 
	 * public MathMap<RationalFunction<E>, SimpleRationalFunction<B>>
	 * getIsomorphism() { return isomorphism; }
	 * 
	 * public MathMap<SimpleRationalFunction<B>, RationalFunction<E>>
	 * getInverseIsomorphism() { return inverseIsomorphism; }
	 * 
	 * }
	 * 
	 * public static <B extends Element<B>, E extends AlgebraicExtensionElement<B,
	 * E>, Ext extends FieldExtension<B, E, Ext>>
	 * SimpleFunctionFieldFromCoordinateRingOverExtension<B, E>
	 * fromCoordinateRingOverExtensionField( Ext extension, CoordinateRing<E>
	 * coordinateRing) { ExtensionCoordinateRing<B, E> extensionCoordinateRing =
	 * extension.asCoordinateRing(coordinateRing);
	 * SimpleFunctionFieldFromCoordinateRing<B> fromBaseRing = fromCoordinateRing(
	 * extensionCoordinateRing.getBaseCoordinateRing()); GenericProjectiveScheme<B>
	 * domainOverBase = GenericProjectiveScheme
	 * .fromAffineCoordinateRing(extensionCoordinateRing.getBaseCoordinateRing());
	 * FunctionField<B> functionField = new FunctionField<>(domainOverBase); return
	 * new SimpleFunctionFieldFromCoordinateRingOverExtension<>(fromBaseRing.
	 * getSimpleFunctionField(), new MathMap<>() {
	 * 
	 * @Override public SimpleRationalFunction<B> evaluate(RationalFunction<E> t) {
	 * Polynomial<B> numerator = extensionCoordinateRing.getIsomorphism()
	 * .evaluate(coordinateRing.getEmbedding(t.getNumerator())).getElement();
	 * Polynomial<B> denominator = extensionCoordinateRing.getIsomorphism()
	 * .evaluate(coordinateRing.getEmbedding(t.getDenominator())).getElement();
	 * return fromBaseRing.isomorphism.evaluate(functionField.getFunction(numerator,
	 * denominator)); } }, new MathMap<>() {
	 * 
	 * @Override public RationalFunction<E> evaluate(SimpleRationalFunction<B> t) {
	 * // Polynomial<B> numerator = extensionCoordinateRing.getIsomorphism() //
	 * .evaluate(coordinateRing.getEmbedding(t.getNumerator())).getElement(); //
	 * Polynomial<B> denominator = extensionCoordinateRing.getIsomorphism() //
	 * .evaluate(coordinateRing.getEmbedding(t.getDenominator())).getElement(); //
	 * return fromBaseRing.isomorphism.evaluate(functionField.getFunction(numerator,
	 * denominator)); // TODO Auto-generated method stub return null; } }); }
	 */

	public static class SimpleFunctionFieldFromCoordinateRing<T extends Element<T>> {
		private SimpleFunctionField<T> simpleFunctionField;
		private MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> isomorphism;
		private MathMap<SimpleRationalFunction<T>, TExt<T>> inverseIsomorphism;

		private SimpleFunctionFieldFromCoordinateRing(SimpleFunctionField<T> simpleFunctionField,
				MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> isomorphism,
				MathMap<SimpleRationalFunction<T>, TExt<T>> inverseIsomorphism) {
			this.simpleFunctionField = simpleFunctionField;
			this.isomorphism = isomorphism;
			this.inverseIsomorphism = inverseIsomorphism;
		}

		public SimpleFunctionField<T> getSimpleFunctionField() {
			return simpleFunctionField;
		}

		public MathMap<CoordinateRingElement<T>, SimpleRationalFunction<T>> getIsomorphism() {
			return isomorphism;
		}

		public MathMap<SimpleRationalFunction<T>, TExt<T>> getInverseIsomorphism() {
			return inverseIsomorphism;
		}
	}

	public static <T extends Element<T>> SimpleFunctionFieldFromCoordinateRing<T> fromCoordinateRing(
			CoordinateRing<T> coordinateRing) {
//		if (!coordinateRing.isIntegral()) {
//			throw new ArithmeticException("Not an integral coordinate ring!");
//		}
		if (!(coordinateRing.getPolynomialRing().getRing() instanceof Field<?>)) {
			throw new ArithmeticException("Not a coordinate ring over a field!");
		}
		Field<T> field = (Field<T>) coordinateRing.getPolynomialRing().getRing();
		String[] freeVariableNames = new String[coordinateRing.krullDimension()];
		String[] boundVariableNames = new String[coordinateRing.boundVariables().size()];
		int freeIndex = 0;
		int boundIndex = 0;
		int[] freeToGeneralMap = new int[coordinateRing.krullDimension()];
		int[] boundToGeneralMap = new int[coordinateRing.boundVariables().size()];
		for (int i = 0; i < coordinateRing.getPolynomialRing().numberOfVariables(); i++) {
			if (coordinateRing.boundVariables().contains(i + 1)) {
				boundVariableNames[boundIndex] = coordinateRing.getPolynomialRing().getVariableNames()[i];
				boundToGeneralMap[boundIndex] = i;
				boundIndex++;
			} else {
				freeVariableNames[freeIndex] = coordinateRing.getPolynomialRing().getVariableNames()[i];
				freeToGeneralMap[freeIndex] = i;
				freeIndex++;
			}
		}
		TranscendentalFieldExtension<T> transcendental = new TranscendentalFieldExtension<>(field, freeVariableNames);
		PolynomialRing<T> basePolynomialRing = transcendental.polynomialRing();
		PolynomialRing<TExt<T>> polynomialRing = AbstractPolynomialRing.getPolynomialRing(transcendental,
				coordinateRing.getPolynomialRing().getComparator(), boundVariableNames);
		MathMap<Polynomial<T>, Polynomial<TExt<T>>> toTranscendentalPolynomial = new FunctionMathMap<>(
				(Polynomial<T> t) -> {
					Map<Monomial, TExt<T>> result = new TreeMap<>();
					for (Monomial m : t.monomials()) {
						int[] exponents = new int[polynomialRing.numberOfVariables()];
						int[] baseExponents = new int[basePolynomialRing.numberOfVariables()];
						int j = 0;
						int k = 0;
						for (int i = 0; i < m.exponents().length; i++) {
							if (coordinateRing.boundVariables().contains(i + 1)) {
								exponents[j] = m.exponents()[i];
								j++;
							} else {
								baseExponents[k] = m.exponents()[i];
								k++;
							}
						}
						Monomial mext = polynomialRing.getMonomial(exponents);
						Monomial mbase = basePolynomialRing.getMonomial(baseExponents);
						if (!result.containsKey(mext)) {
							result.put(mext, transcendental.zero());
						}
						result.put(mext, transcendental.add(result.get(mext), transcendental.getEmbedding(
								basePolynomialRing.getPolynomial(Collections.singletonMap(mbase, t.coefficient(m))))));
					}
					return polynomialRing.getPolynomial(result);
				});
		List<Polynomial<TExt<T>>> overTranscendental = new ArrayList<>();
		for (Polynomial<T> generator : coordinateRing.getIdeal().generators()) {
			overTranscendental.add(toTranscendentalPolynomial.evaluate(generator));
		}
		PolynomialIdeal<TExt<T>> ideal = polynomialRing.getIdeal(overTranscendental);
		int[] degrees = new int[polynomialRing.numberOfVariables()];
		// for (Polynomial<TExt<T>> generator : ideal.generators()) {
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			// int degreeGenerator = generator.degree(i + 1);
			// if (degrees[i] <= degreeGenerator) {
			degrees[i] = ideal.degree() /* degreeGenerator */ + 1;
			// }
		}
		// }
		CoordinateRing<TExt<T>> transcendentalCoordinateRing = ideal.divideOut();
		if (transcendentalCoordinateRing.krullDimension() != 0) {
			throw new ArithmeticException("Algorithm wrong");
		}
		int degree = transcendentalCoordinateRing.degree();
		Integers z = Integers.z();
		UnivariatePolynomial<TExt<T>> minimalPolynomial = null;
		Polynomial<TExt<T>> candidate = null;
		MathMap<CoordinateRingElement<TExt<T>>, Vector<TExt<T>>> asVector = new FunctionMathMap<>(
				(CoordinateRingElement<TExt<T>> t) -> polynomialRing.asVector(t.getElement(), degrees));
		if (field.characteristic().equals(BigInteger.ZERO)) {
			// Ugly hack!
			Iterable<UnivariatePolynomial<IntE>> intIterator = z.getUnivariatePolynomialRing()
					.polynomialSet(polynomialRing.numberOfVariables() - 1);
			for (Polynomial<IntE> it : intIterator) {
				Vector<IntE> vector = z.getUnivariatePolynomialRing().asVector(it,
						polynomialRing.numberOfVariables() - 1);
				candidate = polynomialRing.zero();
				int i = 0;
				for (IntE coefficient : vector.asList()) {
					candidate = polynomialRing.add(candidate,
							polynomialRing.multiply(coefficient, polynomialRing.getVar(i + 1)));
					i++;
				}
				minimalPolynomial = AbstractFieldExtension.minimalPolynomial(
						transcendentalCoordinateRing.getEmbedding(candidate), degree, transcendentalCoordinateRing,
						transcendental, asVector);
				if (minimalPolynomial.degree() == degree) {
					break;
				}
			}
		} else {
			FiniteVectorSpace<T> module = new FiniteVectorSpace<>(field, polynomialRing.numberOfVariables());
			for (Vector<T> vector : module) {
				candidate = polynomialRing.zero();
				int i = 0;
				for (T coefficient : vector.asList()) {
					candidate = polynomialRing.add(candidate, polynomialRing
							.multiply(transcendental.getEmbedding(coefficient), polynomialRing.getVar(i + 1)));
					i++;
				}
				minimalPolynomial = AbstractFieldExtension.minimalPolynomial(
						transcendentalCoordinateRing.getEmbedding(candidate), degree, transcendentalCoordinateRing,
						transcendental, asVector);
				if (minimalPolynomial.degree() == degree) {
					break;
				}
			}
		}
		List<Vector<TExt<T>>> powers = new ArrayList<>();
		for (int i = 0; i < minimalPolynomial.degree(); i++) {
			powers.add(asVector.evaluate(
					transcendentalCoordinateRing.power(transcendentalCoordinateRing.getEmbedding(candidate), i)));
		}
		Matrix<TExt<T>> baseChange = Matrix.fromColumns(powers);
		SimpleFunctionField<T> functionField = new SimpleFunctionField<>(minimalPolynomial, transcendental,
				polynomialRing.getVariableNames()[0]);
		TranscendentalFieldExtension<T> allVariables = new TranscendentalFieldExtension<>(field,
				coordinateRing.getPolynomialRing());
		final Polynomial<TExt<T>> candidateFinal = candidate;
		return new SimpleFunctionFieldFromCoordinateRing<>(functionField, new MathMap<>() {

			@Override
			public SimpleRationalFunction<T> evaluate(CoordinateRingElement<T> t) {
				Polynomial<TExt<T>> overTranscendental = toTranscendentalPolynomial.evaluate(t.getElement());
				Vector<TExt<T>> vector = asVector
						.evaluate(transcendentalCoordinateRing.getEmbedding(overTranscendental));
				Vector<TExt<T>> inPowerBasis = baseChange.getModule(transcendental).solve(baseChange, vector);
				return functionField.fromVector(inPowerBasis);
			}
		}, new MathMap<>() {
			@Override
			public TExt<T> evaluate(SimpleRationalFunction<T> t) {
				CoordinateRingElement<TExt<T>> overTranscendental = transcendentalCoordinateRing.zero();
				CoordinateRingElement<TExt<T>> candidatePower = transcendentalCoordinateRing.one();
				for (int i = 0; i < degree; i++) {
					overTranscendental = transcendentalCoordinateRing.add(
							transcendentalCoordinateRing.multiply(candidatePower,
									transcendentalCoordinateRing
											.getEmbedding(t.asPolynomial().univariateCoefficient(i))),
							overTranscendental);
					candidatePower = transcendentalCoordinateRing
							.multiply(transcendentalCoordinateRing.getEmbedding(candidateFinal), candidatePower);
				}
				Polynomial<T> freeDenominator = transcendental.polynomialRing().one();
				for (Monomial m : overTranscendental.getElement().monomials()) {
					freeDenominator = transcendental.polynomialRing()
							.lcm(overTranscendental.getElement().coefficient(m).getDenominator(), freeDenominator);
				}
				Polynomial<T> denominator = coordinateRing.getPolynomialRing().getEmbedding(freeDenominator,
						freeToGeneralMap);
				Polynomial<TExt<T>> denominatorRemoved = polynomialRing
						.multiply(transcendental.getEmbedding(freeDenominator), overTranscendental.getElement());
				Polynomial<T> numerator = coordinateRing.getPolynomialRing().zero();
				for (Monomial m : denominatorRemoved.monomials()) {
					int[] exponents = new int[coordinateRing.getPolynomialRing().numberOfVariables()];
					for (int i = 0; i < coordinateRing.boundVariables().size(); i++) {
						exponents[boundToGeneralMap[i]] = m.exponents()[i];
					}
					Polynomial<T> monomial = coordinateRing.getPolynomialRing().getEmbedding(field.one(), exponents);
					numerator = coordinateRing.getPolynomialRing().add(
							coordinateRing.getPolynomialRing().multiply(monomial,
									coordinateRing.getPolynomialRing().getEmbedding(
											denominatorRemoved.coefficient(m).asInteger(), freeToGeneralMap)),
							numerator);
				}
				numerator = coordinateRing.getEmbedding(numerator).getElement();
				denominator = coordinateRing.getEmbedding(denominator).getElement();
				return allVariables.getElement(numerator, denominator);
			}
		});
	}

	public static <T extends Element<T>> SimpleFunctionField<T> getSimpleFunctionField(Field<T> field,
			PolynomialRing<T> polynomialRing, Polynomial<T> polynomial, int variable) {
		if (polynomial.degree(variable) == 0) {
			throw new ArithmeticException("Polynomial is trivial in variable " + variable);
		}
		TranscendentalFieldExtension<T> transcendental = new TranscendentalFieldExtension<>(field,
				polynomialRing.eliminateVariable());
		UnivariatePolynomial<Polynomial<T>> asUnivariate = polynomialRing.asUnivariatePolynomial(polynomial, variable);
		UnivariatePolynomial<TExt<T>> minimalPolynomial = transcendental.getUnivariatePolynomialRing()
				.getEmbedding(asUnivariate, new MathMap<>() {
					@Override
					public TExt<T> evaluate(Polynomial<T> t) {
						return transcendental.getEmbedding(t);
					}
				});
		return new SimpleFunctionField<>(minimalPolynomial, transcendental,
				polynomialRing.getVariableNames()[variable - 1]);
	}

	public SimpleFunctionField(UnivariatePolynomial<TExt<T>> minimalPolynomial, TranscendentalFieldExtension<T> base,
			String variableName) {
		super(minimalPolynomial, base, variableName);
		this.transcendentalExtension = base;
	}

	@Override
	protected SimpleRationalFunction<T> fromSmallDegreePolynomial(UnivariatePolynomial<TExt<T>> polynomial) {
		return new SimpleRationalFunction<>(this, polynomial);
	}

	@Override
	public SimpleFunctionField<T> makeExtension(UnivariatePolynomial<TExt<T>> minimalPolynomial) {
		return new SimpleFunctionField<>(minimalPolynomial, transcendentalExtension, getVariableName());
	}

	@Override
	protected SimpleFunctionField<T> asExtensionType() {
		return this;
	}

	@Override
	public TranscendentalFieldExtension<T> getBaseField() {
		return transcendentalExtension;
	}

	@Override
	public boolean hasCharacteristicRoot(SimpleRationalFunction<T> t, int power) {
		return getBaseField().getUnivariatePolynomialRing().hasCharacteristicRoot(t.asPolynomial());
	}

	@Override
	public SimpleRationalFunction<T> characteristicRoot(SimpleRationalFunction<T> t, int power) {
		UnivariatePolynomialRing<TExt<T>> polynomialRing = getBaseField().getUnivariatePolynomialRing();
		return fromPolynomial(polynomialRing.toUnivariate(polynomialRing.characteristicRoot(t.asPolynomial())));
	}
}
