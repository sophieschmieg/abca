package varieties.projective;

import java.util.ArrayList;
import java.util.List;

import fields.helper.FieldOfFractions.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.local.Value;
import fields.local.ValueGroup;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.polynomials.PolynomialIdeal;
import varieties.AbstractMorphism;
import varieties.FunctionField;
import varieties.GeneralRationalFunction;
import varieties.RationalFunction;
import varieties.affine.AffineMorphism;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;

public class ProjectiveMorphism<T extends Element<T>>
		extends AbstractMorphism<T, ProjectivePoint<T>, ProjectivePoint<T>> {
	private Field<T> field;
	private GenericProjectiveScheme<T> domain;
	private GenericProjectiveScheme<T> range;
	private List<Polynomial<T>> asPolynomials;
	private PolynomialRing<T> polynomials;

	public ProjectiveMorphism(GenericProjectiveScheme<T> domain, GenericProjectiveScheme<T> range,
			List<Polynomial<T>> asPolynomials) {
		if (!domain.getField().equals(range.getField())) {
			throw new ArithmeticException("Fields mismatched!");
		}
		this.field = domain.getField();
		this.domain = domain;
		this.range = range;
		int degree = -1;
		this.polynomials = domain.homogenousPolynomialRing();
		Polynomial<T> gcd = polynomials.zero();
		for (Polynomial<T> polynomial : asPolynomials) {
			if (polynomial.numberOfVariables() != polynomials.numberOfVariables()) {
				throw new ArithmeticException("domain dimensions mismatched!");
			}
			if (degree == -1) {
				degree = polynomial.degree();
			}
			if (!polynomials.isHomogeneous(polynomial)) {
				throw new ArithmeticException("polynomials not homogenous");
			}
			if (degree != polynomial.degree() && !polynomial.equals(polynomials.zero())) {
				throw new ArithmeticException("polynomials not of equal degree");
			}
			gcd = polynomials.gcd(polynomial, gcd);
		}
		if (degree == -1) {
			throw new ArithmeticException("No nonzero polynomial!");
		}
		if (asPolynomials.size() != range.homogenousPolynomialRing().numberOfVariables()) {
			throw new ArithmeticException("range dimensions mismatched!");
		}
		this.asPolynomials = new ArrayList<>();
		for (Polynomial<T> polynomial : asPolynomials) {
			this.asPolynomials.add(polynomials.divideChecked(polynomial, gcd));
		}
	}

	public static <T extends Element<T>> ProjectiveMorphism<T> fromRationalFunctions(GenericProjectiveScheme<T> domain,
			GenericProjectiveScheme<T> range, List<RationalFunction<T>> asRationalFunctions, int rangeAffineIndex) {
		FunctionField<T> functionField = domain.getFunctionField();
		CoordinateRing<T> cr = functionField.affineCoordinateRing();
		PolynomialRing<T> polynomialRing = cr.getPolynomialRing();
		Polynomial<T> lcm = polynomialRing.one();
		for (RationalFunction<T> rationalFunction : asRationalFunctions) {
			lcm = polynomialRing.lcm(functionField.getIntegerDenominator(rationalFunction).getElement(), lcm);
		}
		lcm = cr.getEmbedding(lcm).getElement();
		int degree = lcm.degree();
		List<Polynomial<T>> reduced = new ArrayList<>();
		for (RationalFunction<T> rationalFunction : asRationalFunctions) {
			Polynomial<T> p = cr.divideChecked(
					cr.multiply(cr.getEmbedding(lcm), functionField.getIntegerNumerator(rationalFunction)),
					functionField.getIntegerDenominator(rationalFunction)).getElement();
			reduced.add(p);
			degree = Math.max(p.degree(), degree);
		}
		List<Polynomial<T>> asPolynomials = new ArrayList<>();
		int index = 0;
		for (int i = 0; i < range.homogenousPolynomialRing().numberOfVariables(); i++) {
			Polynomial<T> p;
			if (i == rangeAffineIndex) {
				p = lcm;
			} else {
				p = reduced.get(index);
				index++;
			}
			Polynomial<T> homogenous = domain.homogenize(p, functionField.affineCoverIndex());
			homogenous = domain.asGenericProjectiveScheme().homogenousPolynomialRing()
					.multiply(domain.asGenericProjectiveScheme().homogenousPolynomialRing().getVarPower(
							functionField.affineCoverIndex() + 1, degree - homogenous.degree()), homogenous);
			asPolynomials.add(homogenous);
		}
		return new ProjectiveMorphism<>(domain, range, asPolynomials);
	}

	public List<Polynomial<T>> asPolynomials() {
		return asPolynomials;
	}

	public boolean isFinite() {
		return false;
	}

	public boolean isInjective() {
		return false;
	}

	private class UpToUniformizerResult {
		private List<LocalizedElement<T>> upToUniformizer;
		private List<Value> values;
		private Value minimumValue;
	}

	private UpToUniformizerResult upToUniformizer(int coverIndex, PolynomialRing<T> affineRing,
			LocalizedCoordinateRing<T> localizedRing) {
		int[] map = new int[domain.homogenousPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map.length; i++) {
			if (i == coverIndex) {
				map[i] = -1;
				continue;
			}
			map[i] = i > coverIndex ? i - 1 : i;
		}
		UpToUniformizerResult result = new UpToUniformizerResult();
		result.upToUniformizer = new ArrayList<>();
		result.values = new ArrayList<>();
		result.minimumValue = Value.INFINITY;
		for (Polynomial<T> polynomial : asPolynomials) {
			Polynomial<T> affinePolynomial = affineRing
					.getEmbedding(domain.homogenousPolynomialRing().dehomogenize(polynomial, coverIndex + 1), map);
			LocalizedElement<T> inLocalization = localizedRing.getEmbedding(affinePolynomial);
			Value value = localizedRing.valuation(inLocalization);
			if (result.minimumValue.compareTo(value) > 0) {
				result.minimumValue = value;
			}
			result.values.add(value);
			result.upToUniformizer.add(localizedRing.upToUniformizerPower(inLocalization));
		}
		return result;
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> t) {
		int coverIndex = domain.affineCoverIndex(t).get(0);
		AffineScheme<T> affineDomain = domain.getAffineCover().getCover().get(coverIndex);
		PolynomialRing<T> affineRing = affineDomain.getCoordinateRing().getPolynomialRing();
		AffinePoint<T> affinePreimage = domain.asAffinePoint(t, coverIndex);
		LocalizedCoordinateRing<T> localizedRing = affineDomain.localizedCoordinateRing(affinePreimage);
		UpToUniformizerResult result = upToUniformizer(coverIndex, affineRing, localizedRing);
		List<T> coords = new ArrayList<>();
		for (int i = 0; i < result.values.size(); i++) {
			if (result.values.get(i).equals(result.minimumValue)) {
				coords.add(localizedRing.reduceInteger(result.upToUniformizer.get(i)));
			} else {
				coords.add(field.zero());
			}
		}
		return new ProjectivePoint<>(field, coords);
	}

	@Override
	public GenericProjectiveScheme<T> getDomain() {
		return domain;
	}

	@Override
	public GenericProjectiveScheme<T> getRange() {
		return range;
	}

	@Override
	public GeneralRationalFunction<T, ProjectivePoint<T>, ProjectivePoint<T>> restrict(ProjectivePoint<T> preimage) {
		ProjectivePoint<T> image = evaluate(preimage);
		int domainCoverIndex = domain.affineCoverIndex(preimage).get(0);
		int rangeCoverIndex = range.affineCoverIndex(image).get(0);
		AffineScheme<T> affineDomain = domain.getAffineCover().getCover().get(domainCoverIndex);
		PolynomialRing<T> affineRing = affineDomain.getCoordinateRing().getPolynomialRing();
		AffinePoint<T> affinePreimage = domain.asAffinePoint(preimage, domainCoverIndex);
		LocalizedCoordinateRing<T> localizedRing = affineDomain.localizedCoordinateRing(affinePreimage);
		UpToUniformizerResult upToUniformizer = upToUniformizer(domainCoverIndex, affineRing, localizedRing);
		if (!upToUniformizer.minimumValue.equals(upToUniformizer.values.get(rangeCoverIndex))) {
			throw new ArithmeticException("Range Cover Index wrong!");
		}
		Fraction<Polynomial<T>> rangeCoverFraction = upToUniformizer.upToUniformizer.get(rangeCoverIndex)
				.asPolynomialFraction();
		rangeCoverFraction.canonicalize();
		LocalizedElement<T> localDenominator = localizedRing.getEmbedding(rangeCoverFraction.getDenominator());

		Polynomial<T> denominator = rangeCoverFraction.getNumerator();
		List<Polynomial<T>> numerators = new ArrayList<>();
		List<Polynomial<T>> denominators = new ArrayList<>();
		ValueGroup g = ValueGroup.g();
		for (int i = 0; i < upToUniformizer.values.size(); i++) {
			if (i == rangeCoverIndex) {
				continue;
			}
			if (upToUniformizer.values.get(i).isInfinite()) {
				numerators.add(affineRing.zero());
				continue;
			}
			LocalizedElement<T> reduced = localizedRing.multiply(localDenominator,
					upToUniformizer.upToUniformizer.get(i), localizedRing.power(localizedRing.uniformizer(),
							g.subtract(upToUniformizer.values.get(i), upToUniformizer.minimumValue).value()));
			Fraction<Polynomial<T>> asFraction = reduced.asPolynomialFraction();
			asFraction.canonicalize();
			denominator = affineRing.lcm(asFraction.getDenominator(), denominator);
			numerators.add(asFraction.getNumerator());
			denominators.add(asFraction.getDenominator());
		}
		Polynomial<T> inverse = polynomials.subtract(polynomials.multiply(polynomials.getEmbedding(denominator),
				polynomials.getVar(polynomials.numberOfVariables())), polynomials.one());
		List<Polynomial<T>> inverseList = new ArrayList<>();
		inverseList.add(inverse);
		for (Polynomial<T> generator : affineDomain.getCoordinateRing().getIdeal().generators()) {
			inverseList.add(polynomials.getEmbedding(generator));
		}
		PolynomialIdeal<T> inverseIdeal = polynomials.getIdeal(inverseList);
		CoordinateRing<T> domainSubsetRing = inverseIdeal.divideOut();
		AffineScheme<T> domainSubset = new AffineScheme<>(field, domainSubsetRing);
		List<CoordinateRingElement<T>> embeddingPolynomials = new ArrayList<>();
		for (int i = 0; i < affineRing.numberOfVariables(); i++) {
			embeddingPolynomials.add(domainSubsetRing.getEmbedding(polynomials.getVar(i + 1)));
		}
		AffineMorphism<T> embedding = new AffineMorphism<>(domainSubset, affineDomain, embeddingPolynomials);
		List<CoordinateRingElement<T>> asPolynomials = new ArrayList<>();
		for (int i = 0; i < numerators.size(); i++) {
			Polynomial<T> extend = affineRing.divideChecked(denominator, denominators.get(i));
			Polynomial<T> numerator = affineRing.multiply(numerators.get(i), extend);
			Polynomial<T> polynomial = polynomials.multiply(polynomials.getEmbedding(numerator),
					polynomials.getVar(polynomials.numberOfVariables()));
			asPolynomials.add(domainSubsetRing.getEmbedding(polynomial));
		}
		AffineMorphism<T> restricted = new AffineMorphism<>(domainSubset, affineDomain, asPolynomials);
		return new GeneralRationalFunction<>(domain, range, restricted, domainCoverIndex, embedding, rangeCoverIndex);
	}
}
