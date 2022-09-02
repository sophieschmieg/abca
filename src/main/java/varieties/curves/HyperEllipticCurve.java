package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import util.MiscAlgorithms;
import varieties.Morphism;
import varieties.RationalFunction;
import varieties.affine.AffinePoint;
import varieties.projective.AbstractProjectiveScheme;
import varieties.projective.GenericProjectiveScheme;
import varieties.projective.ProjectivePoint;

public class HyperEllipticCurve<T extends Element<T>> extends AbstractProjectiveScheme<T>
		implements SmoothCurve<T>, Element<HyperEllipticCurve<T>> {
	private Field<T> field;
	private Polynomial<T> definingPolynomial;
	private PolynomialRing<T> affinePolynomialRing;
	private UnivariatePolynomial<T> rhs;
	private UnivariatePolynomial<T> yFactor;
	private int genus;
	private GenericProjectiveScheme<T> asGenericProjectiveScheme;

	public HyperEllipticCurve(Field<T> field, UnivariatePolynomial<T> rhs, UnivariatePolynomial<T> yFactor) {
		this.field = field;
		UnivariatePolynomialRing<T> univariatePolynomials = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<T> rhsDerivative = univariatePolynomials.derivative(rhs);
		UnivariatePolynomial<T> yFactorDerivative = univariatePolynomials.derivative(yFactor);
		if (field.characteristic().equals(BigInteger.TWO)) {
			if (yFactor.degree() > 0 || yFactor.equals(univariatePolynomials.zero())) {
				UnivariatePolynomial<T> adjusted = univariatePolynomials.add(
						univariatePolynomials.multiply(rhsDerivative, rhsDerivative),
						univariatePolynomials.multiply(-1, rhs, yFactorDerivative, yFactorDerivative));
				if (univariatePolynomials.gcd(yFactor, adjusted).degree() != 0) {
					throw new ArithmeticException("Singular curve");
				}
			}
			List<T> lcRelationList = new ArrayList<>();
			lcRelationList.add(rhs.leadingCoefficient());
			lcRelationList.add(yFactor.leadingCoefficient());
			lcRelationList.add(field.one());
			UnivariatePolynomial<T> lcRelation = univariatePolynomials.getPolynomial(lcRelationList);
			while (rhs.degree() > 2 * yFactor.degree() && rhs.degree() % 2 == 0 || yFactor.degree() <= rhs.degree()
					|| (rhs.degree() == 2 * yFactor.degree() && field.hasRoots(lcRelation))) {
				T coeff;
				int degree;
				if (rhs.degree() > 2 * yFactor.degree() && rhs.degree() % 2 == 0) {
					coeff = field.characteristicRoot(rhs.leadingCoefficient());
					degree = rhs.degree() / 2;
				} else if (yFactor.degree() <= rhs.degree()) {
					coeff = field.divide(rhs.leadingCoefficient(), yFactor.leadingCoefficient());
					degree = rhs.degree() - yFactor.degree();
				} else if (rhs.degree() == 2 * yFactor.degree() && field.hasRoots(lcRelation)) {
					coeff = field.roots(lcRelation).keySet().iterator().next();
					degree = rhs.degree() / 2;
				} else {
					throw new ArithmeticException("If statements wrong!");
				}
				UnivariatePolynomial<T> yOffset = univariatePolynomials.getEmbedding(coeff, degree);
				rhs = univariatePolynomials
						.toUnivariate(univariatePolynomials.add(rhs, univariatePolynomials.multiply(yOffset, yOffset),
								univariatePolynomials.multiply(yOffset, yFactor)));
				lcRelationList.clear();
				lcRelationList.add(rhs.leadingCoefficient());
				lcRelationList.add(yFactor.leadingCoefficient());
				lcRelationList.add(field.one());
				lcRelation = univariatePolynomials.getPolynomial(lcRelationList);
			}
		} else {
			UnivariatePolynomial<T> adjusted = univariatePolynomials
					.add(univariatePolynomials.multiply(yFactor, yFactor), univariatePolynomials.multiply(4, rhs));
			UnivariatePolynomial<T> adjustedDerivative = univariatePolynomials.add(
					univariatePolynomials.multiply(yFactor, yFactorDerivative),
					univariatePolynomials.multiply(2, rhsDerivative));
			if (univariatePolynomials.gcd(adjusted, adjustedDerivative).degree() != 0) {
				throw new ArithmeticException("Singular curve");
			}
			int degree = Math.max(MiscAlgorithms.DivRoundUp(rhs.degree(), 2), yFactor.degree());
			T half = field.inverse(field.getInteger(2));
			T yFactorDegree = yFactor.univariateCoefficient(degree);
			T yFactorDegreeMinus = yFactor.univariateCoefficient(degree - 1);
			T rhs2Degree = rhs.univariateCoefficient(2 * degree);
			T rhs2DegreeMinus = rhs.univariateCoefficient(2 * degree - 1);
			while (!yFactorDegree.equals(field.zero())
					&& rhs2Degree.equals(field.negative(field.power(field.multiply(half, yFactorDegree), 2)))
					&& rhs2DegreeMinus.equals(field.multiply(-1, half, yFactorDegree, yFactorDegreeMinus))) {
				UnivariatePolynomial<T> yOffset = univariatePolynomials
						.getEmbedding(field.multiply(-1, half, yFactorDegree), degree);
				UnivariatePolynomial<T> yAdjustment = univariatePolynomials
						.toUnivariate(univariatePolynomials.multiply(2, yOffset));
				UnivariatePolynomial<T> rhsAdjustment = univariatePolynomials
						.negative(univariatePolynomials.add(univariatePolynomials.multiply(yOffset, yOffset),
								univariatePolynomials.multiply(yFactor, yOffset)));
				rhs = univariatePolynomials.add(rhs, rhsAdjustment);
				yFactor = univariatePolynomials.add(yFactor, yAdjustment);
				degree = Math.max(MiscAlgorithms.DivRoundUp(rhs.degree(), 2), yFactor.degree());
				yFactorDegree = yFactor.univariateCoefficient(degree);
				yFactorDegreeMinus = yFactor.univariateCoefficient(degree - 1);
				rhs2Degree = rhs.univariateCoefficient(2 * degree);
				rhs2DegreeMinus = rhs.univariateCoefficient(2 * degree - 1);
			}
		}
		this.rhs = rhs;
		this.yFactor = yFactor;
		this.genus = Math.max(yFactor.degree() - 1, MiscAlgorithms.DivRoundDown(rhs.degree() - 1, 2));
		if (genus < 2) {
			throw new ArithmeticException("Not a hyperelliptic curve!");
		}
		affinePolynomialRing = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.GREVLEX);
		this.definingPolynomial = affinePolynomialRing.add(affinePolynomialRing.getVarPower(2, 2),
				affinePolynomialRing.multiply(affinePolynomialRing.getEmbedding(yFactor),
						affinePolynomialRing.getVar(2)),
				affinePolynomialRing.negative(affinePolynomialRing.getEmbedding(rhs)));

	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		if (p.getDim() != genus + 2) {
			return false;
		}
		if (!p.getCoord(2).equals(field.zero())) {
			for (int i = 0; i < genus + 2; i++) {
				if (!p.getDehomogenisedCoord(i + 1, 2).equals(field.one())) {
					return false;
				}
			}
			AffinePoint<T> affine = new AffinePoint<>(field, p.getDehomogenisedCoord(2, 1),
					p.getDehomogenisedCoord(genus + 3, 1));
			return affinePolynomialRing.evaluate(definingPolynomial, affine).equals(field.zero());
		}
		for (int i = 0; i < genus; i++) {
			if (!p.getCoord(i + 1).equals(field.zero())) {
				return false;
			}
		}
		T y = p.getCoord(genus + 3);
		T x = p.getCoord(genus + 2);
		return field.add(field.multiply(y, y), field.multiply(yFactor.univariateCoefficient(genus + 1), y, x),
				field.multiply(-1, rhs.univariateCoefficient(2 * genus + 2), x, x)).equals(field.zero());
	}

	@Override
	public Field<T> getField() {
		return this.field;
	}

	@Override
	public int getEmbeddingDimension() {
		return genus + 2;
	}

	@Override
	public String toString() {
		return this.definingPolynomial.toString();
	}

	@Override
	public Morphism<T, ProjectivePoint<T>, ProjectivePoint<T>> identityMorphism() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<? extends Morphism<T, ProjectivePoint<T>, ProjectivePoint<T>>> irreducibleComponents() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Morphism<T, ProjectivePoint<T>, ProjectivePoint<T>> reduced() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFinite() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterator<ProjectivePoint<T>> iterator() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int compareTo(HyperEllipticCurve<T> o) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public GenericProjectiveScheme<T> asGenericProjectiveScheme() {
		if (asGenericProjectiveScheme == null) {
			PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, genus + 3, Monomial.GREVLEX);
			List<Polynomial<T>> equations = new ArrayList<>();
			Polynomial<T> defining = r.getVarPower(genus + 3, 2);
			for (int i = 0; i <= genus + 1; i++) {
				defining = r.add(r.multiply(yFactor.univariateCoefficient(i), r.getVar(i + 1), r.getVar(genus + 3)),
						defining);
				defining = r.add(r.multiply(rhs.univariateCoefficient(2 * i), r.getVar(i + 1), r.getVar(i + 1)),
						defining);
				defining = r.add(r.multiply(rhs.univariateCoefficient(2 * i + 1), r.getVar(i + 1), r.getVar(i + 2)),
						defining);
			}
			equations.add(defining);
			for (int i = 1; i < genus + 1; i++) {
				equations.add(r.subtract(r.getVarPower(i + 1, 2), r.multiply(r.getVar(i), r.getVar(i + 2))));
				equations.add(r.subtract(r.multiply(r.getVar(i), r.getVar(genus + 2)),
						r.multiply(r.getVar(MiscAlgorithms.DivRoundUp(genus + i, 2) + 1),
								r.getVar(MiscAlgorithms.DivRoundDown(genus + i, 2) + 1))));
			}
			asGenericProjectiveScheme = new GenericProjectiveScheme<>(field, r, equations);
		}
		return asGenericProjectiveScheme;
	}

	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(WeilDivisor<T> div) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isPrincipal(WeilDivisor<T> div) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public CoordinateRing<T> getCoordinateRing() {
		// TODO Auto-generated method stub
		return null;
	}
}
