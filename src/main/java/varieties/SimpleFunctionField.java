package varieties;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.CoordinateRing;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.FreeModule;
import fields.vectors.Vector;
import varieties.SimpleFunctionField.SimpleRationalFunction;
import varieties.affine.AffineScheme;
import varieties.projective.ProjectiveScheme;

public class SimpleFunctionField<T extends Element<T>>
		extends AbstractFieldExtension<TExt<T>, SimpleRationalFunction<T>, SimpleFunctionField<T>> {
	public static class SimpleRationalFunction<T extends Element<T>> extends AbstractElement<SimpleRationalFunction<T>>
			implements AlgebraicExtensionElement<TExt<T>, SimpleRationalFunction<T>> {
		private UnivariatePolynomial<TExt<T>> asPolynomial;

		private SimpleRationalFunction(UnivariatePolynomial<TExt<T>> asPolynomial) {
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
	}

	private TranscendentalFieldExtension<T> transcendentalExtension;
	
	public static <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>> SimpleFunctionField<B> forProjectiveVarietyOverExtensionField(
			Ext extension, ProjectiveScheme<E> domain) {
		int affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
		AffineScheme<E> affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
		return fromCoordinateRingOverExtensionField(extension, affineSlice.getCoordinateRing());
	}

	public static <T extends Element<T>> SimpleFunctionField<T> forProjectiveVariety(
			ProjectiveScheme<T> domain) {
		int affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
		AffineScheme<T> affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
		return fromCoordinateRing(affineSlice.getCoordinateRing());
	}

	public static <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, Ext extends FieldExtension<B, E, Ext>> SimpleFunctionField<B> fromCoordinateRingOverExtensionField(
			Ext extension, CoordinateRing<E> coordinateRing) {
		return fromCoordinateRing(extension.asCoordinateRing(coordinateRing).getBaseCoordinateRing());
	}

	public static <T extends Element<T>> SimpleFunctionField<T> fromCoordinateRing(CoordinateRing<T> coordinateRing) {
		if (!coordinateRing.isIntegral()) {
			throw new ArithmeticException("Not an integral coordinate ring!");
		}
		if (!(coordinateRing.getPolynomialRing().getRing() instanceof Field<?>)) {
			throw new ArithmeticException("Not a coordinate ring over a field!");
		}
		Field<T> field = (Field<T>) coordinateRing.getPolynomialRing().getRing();
		TranscendentalFieldExtension<T> transcendental = new TranscendentalFieldExtension<>(field,
				coordinateRing.krullDimension());
		PolynomialRing<T> basePolynomialRing = transcendental.polynomialRing();
		PolynomialRing<TExt<T>> polynomialRing = AbstractPolynomialRing.getPolynomialRing(transcendental,
				coordinateRing.getPolynomialRing().numberOfVariables() - transcendental.transcendenceDegree(),
				coordinateRing.getPolynomialRing().getComparator());
		MathMap<Polynomial<T>, Polynomial<TExt<T>>> toTranscendentalPolynomial = new MathMap<>() {
			@Override
			public Polynomial<TExt<T>> evaluate(Polynomial<T> t) {
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
			}
		};
		List<Polynomial<TExt<T>>> overTranscendental = new ArrayList<>();
		for (Polynomial<T> generator : coordinateRing.getIdeal().generators()) {
			overTranscendental.add(toTranscendentalPolynomial.evaluate(generator));
		}
		PolynomialIdeal<TExt<T>> ideal = polynomialRing.getIdeal(overTranscendental);
		int[] degrees = new int[polynomialRing.numberOfVariables()];
		for (Polynomial<TExt<T>> generator : ideal.generators()) {
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				int degreeGenerator = generator.degree(i + 1);
				if (degrees[i] < degreeGenerator) {
					degrees[i] = degreeGenerator;
				}
			}
		}
		CoordinateRing<TExt<T>> transcendentalCoordinateRing = new CoordinateRing<>(polynomialRing, ideal);
		if (transcendentalCoordinateRing.krullDimension() != 0) {
			throw new ArithmeticException("Algorithm wrong");
		}
		int degree = transcendentalCoordinateRing.degree();
		Integers z = Integers.z();
		UnivariatePolynomial<TExt<T>> minimalPolynomial = null;
		if (field.characteristic().equals(BigInteger.ZERO)) {
			FreeModule<IntE> intModule = new FreeModule<>(z, polynomialRing.numberOfVariables());
			for (Vector<IntE> vector : intModule) {
				Polynomial<TExt<T>> candidate = polynomialRing.zero();
				int i = 0;
				for (IntE coefficient : vector.asList()) {
					candidate = polynomialRing.add(candidate,
							polynomialRing.multiply(coefficient, polynomialRing.getVar(i + 1)));
					i++;
				}
				minimalPolynomial = AbstractFieldExtension.minimalPolynomial(candidate, degree, polynomialRing,
						transcendental, new MathMap<>() {
							@Override
							public Vector<TExt<T>> evaluate(Polynomial<TExt<T>> t) {
								return polynomialRing.asVector(t, degrees);
							}
						});
				if (minimalPolynomial.degree() == degree) {
					break;
				}
			}
		} else {
			FiniteVectorSpace<T> module = new FiniteVectorSpace<>(field,
					coordinateRing.getPolynomialRing().numberOfVariables());
			for (Vector<T> vector : module) {
				Polynomial<T> candidate = coordinateRing.getPolynomialRing().zero();
				int i = 0;
				for (T coefficient : vector.asList()) {
					candidate = coordinateRing.getPolynomialRing().add(candidate, coordinateRing.getPolynomialRing()
							.multiply(coefficient, coordinateRing.getPolynomialRing().getVar(i + 1)));
					i++;
				}
				minimalPolynomial = AbstractFieldExtension.minimalPolynomial(
						toTranscendentalPolynomial.evaluate(candidate), degree, polynomialRing, transcendental,
						new MathMap<>() {
							@Override
							public Vector<TExt<T>> evaluate(Polynomial<TExt<T>> t) {
								return polynomialRing.asVector(t, degrees);
							}
						});
				if (minimalPolynomial.degree() == degree) {
					break;
				}
			}
		}
		return new SimpleFunctionField<>(minimalPolynomial, transcendental);
	}

	public SimpleFunctionField<T> getSimpleFunctionField(Field<T> field, PolynomialRing<T> polynomialRing,
			Polynomial<T> polynomial, int variable) {
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
		return new SimpleFunctionField<>(minimalPolynomial, transcendental);
	}

	public SimpleFunctionField(UnivariatePolynomial<TExt<T>> minimalPolynomial, TranscendentalFieldExtension<T> base) {
		super(minimalPolynomial, base);
		this.transcendentalExtension = base;
	}

	@Override
	protected SimpleRationalFunction<T> fromSmallDegreePolynomial(UnivariatePolynomial<TExt<T>> polynomial) {
		return new SimpleRationalFunction<>(polynomial);
	}

	@Override
	public SimpleFunctionField<T> makeExtension(UnivariatePolynomial<TExt<T>> minimalPolynomial) {
		return new SimpleFunctionField<>(minimalPolynomial, transcendentalExtension);
	}

	@Override
	protected SimpleFunctionField<T> asExtensionType() {
		return this;
	}
}
