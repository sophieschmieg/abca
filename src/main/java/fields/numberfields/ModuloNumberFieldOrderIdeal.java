package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractRing;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DiscreteValuationField.OtherVersion;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Ring;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.ModuloNumberFieldOrderIdeal.ModNFE;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldOrder.NumberFieldOrderIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.GenericPIDModule;
import fields.vectors.GenericPIDModule.Mod;
import fields.vectors.RationalLattice;
import fields.vectors.SubLattice;
import fields.vectors.Vector;

public class ModuloNumberFieldOrderIdeal extends AbstractRing<ModNFE> implements Ring<ModNFE> {
	public static class ModNFE extends AbstractElement<ModNFE> {
		private NFE lift;

		private ModNFE(NFE lift) {
			this.lift = lift;
		}

		@Override
		public int compareTo(ModNFE o) {
			return lift.compareTo(o.lift);
		}

		@Override
		public String toString() {
			return lift.toString();
		}
	}

	private NumberFieldOrder order;
	private NumberFieldOrderIdeal ideal;
	private RationalLattice orderAsRationalLattice;
	private SubLattice<Vector<Fraction>, Fraction, Vector<Fraction>> idealAsSubLattice;
	private GenericPIDModule<IntE, NFE> genericIntModule;
	private List<ModNFE> unitGenerators;
//	private GenericPIDModule<IntE, Vector<IntE>> genericModule;
//	private Map<ModNFE, Vector<IntE>> knownVectors;
//	private List<ModuloNumberFieldOrderIdeal> primePowerDecomposition;

	ModuloNumberFieldOrderIdeal(NumberFieldOrder order, NumberFieldOrderIdeal ideal) {
		Rationals q = Rationals.q();
		this.order = order;
		this.ideal = ideal;
		List<NFE> reducedIdealBasis = new SubLattice<>(order, ideal.asSubModule(), 1.0).getModuleGenerators();
		this.orderAsRationalLattice = RationalLattice.fromIntegerLattice(order.rank(),
				new FreeModule<>(Integers.z(), order.rank()).getBasis(), 1.0);
		List<Vector<Fraction>> idealGeneratorsIntVectors = new ArrayList<>();
		for (NFE basis : reducedIdealBasis) {
			idealGeneratorsIntVectors.add(Vector.mapVector(q.getEmbeddingMap(), order.asVector(basis)));
		}
		this.idealAsSubLattice = new SubLattice<>(orderAsRationalLattice, idealGeneratorsIntVectors, 1.0);
		List<NFE> reducedBasis = new ArrayList<>();
		for (Vector<Fraction> basisVector : idealAsSubLattice.getModuleGenerators()) {
			Vector<IntE> intVector = Vector.mapVector(q.getAsIntegerMap(), basisVector);
			reducedBasis.add(order.fromVector(intVector));
		}
		this.genericIntModule = new GenericPIDModule<>(order, new FreeSubModule<>(order, reducedBasis));
//		this.knownVectors = new TreeMap<>();
	}

	@Override
	public String toString() {
		return order + "/" + ideal;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public ModNFE getRandomElement() {
		return fromGeneric(genericIntModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return genericIntModule.getNumberOfElements();
	}

	@Override
	public BigInteger getNumberOfUnits() throws InfinityException {
//		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(ideal);
		BigInteger result = BigInteger.ONE;
//		for (Ideal<NFE> primeFactor : factorization.primeFactors()) {
//			NumberFieldIdeal prime = (NumberFieldIdeal) primeFactor;
//			BigInteger q = prime.norm().getValue();
//			result = result.multiply(q.subtract(BigInteger.ONE)).multiply(q.pow(factorization.multiplicity(prime) - 1));
//		}
		return result;
	}

	@Override
	public Iterator<ModNFE> iterator() {
		return new Iterator<>() {
			Iterator<Mod<NFE>> it = genericIntModule.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public ModNFE next() {
				return fromGeneric(it.next());
			}

		};
	}

	@Override
	public ModNFE zero() {
		return reduce(order.zero());
	}

	@Override
	public ModNFE one() {
		return reduce(order.one());
	}

	@Override
	public BigInteger characteristic() {
		return null;// return ideal.intGenerator().getValue();
	}

	public ModNFE reduce(NFE t) {
		Rationals q = Rationals.q();
		NFE reduced = genericIntModule.lift(genericIntModule.reduce(t));
		Vector<Fraction> asFractions = Vector.mapVector(q.getEmbeddingMap(), order.asVector(reduced));
		Vector<Fraction> closestLatticePoint = orderAsRationalLattice.getVectorSpace().closestLatticePoint(asFractions,
				idealAsSubLattice);
		Vector<Fraction> diff = orderAsRationalLattice.getVectorSpace().subtract(asFractions, closestLatticePoint);
		return new ModNFE(order.fromVector(Vector.mapVector(q.getAsIntegerMap(), diff)));
	}

	public NFE lift(ModNFE t) {
		return t.lift;
	}

	private ModNFE fromGeneric(Mod<NFE> t) {
		return reduce(genericIntModule.lift(t));
	}

//	private Mod<NFE> asGeneric(ModNFE t) {
//		return genericIntModule.reduce(lift(t));
//	}

	@Override
	public ModNFE add(ModNFE t1, ModNFE t2) {
		return reduce(order.add(lift(t1), lift(t2)));
	}

	@Override
	public ModNFE negative(ModNFE t) {
		return reduce(order.negative(lift(t)));
	}

	@Override
	public ModNFE multiply(ModNFE t1, ModNFE t2) {
		return reduce(order.multiply(lift(t1), lift(t2)));
	}

	@Override
	public boolean isUnit(ModNFE t) {
		return order.coprime(ideal, order.getIdeal(Collections.singletonList(lift(t))));
	}

	@Override
	public ModNFE inverse(ModNFE t) {
		return reduce(
				order.bezoutIdentity(order.getIdeal(Collections.singletonList(lift(t))), ideal).getCoeff1().get(0));
	}

	@Override
	public boolean isCommutative() {
		return true;
	}

	@Override
	public boolean isIntegral() {
		return ideal.isMaximal();
	}

	@Override
	public boolean isReduced() {
		return order.idealFactorization(ideal).squareFree();
	}

	@Override
	public boolean isIrreducible() {
		return order.idealFactorization(ideal).primeFactors().size() == 1;
	}

	@Override
	public boolean isZeroDivisor(ModNFE t) {
		return !isUnit(t);
	}

	@Override
	public boolean isEuclidean() {
		return isIntegral();
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return isIntegral();
	}

	@Override
	public FactorizationResult<ModNFE, ModNFE> uniqueFactorization(ModNFE t) {
		if (!isIntegral()) {
			throw new ArithmeticException("Not a UFD");
		}
		return new FactorizationResult<>(t, Collections.emptySortedMap());
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return isIntegral();
	}

	@Override
	public boolean isDedekindDomain() {
		return isIntegral();
	}

	// t * i = a + a1 * b1 + a2 * b2 + ... an * bn
	@Override
	public boolean isDivisible(ModNFE dividend, ModNFE divisor) {
		return order.add(order.getIdeal(Collections.singletonList(lift(divisor))), ideal).contains(lift(dividend));
	}

	@Override
	public QuotientAndRemainderResult<ModNFE> quotientAndRemainder(ModNFE dividend, ModNFE divisor) {
		List<NFE> generators = new ArrayList<>();
		generators.add(divisor.lift);
		generators.addAll(ideal.generators());
		IdealResult<NFE, NumberFieldOrderIdeal> asIdeal = order.getIdealWithTransforms(generators);
		List<NFE> generate = asIdeal.getIdeal().generate(dividend.lift);
		ModNFE residue = reduce(asIdeal.getIdeal().residue(dividend.lift));
		ModNFE quotient = zero();
		for (int i = 0; i < generate.size(); i++) {
			quotient = add(quotient,
					reduce(order.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(0))));
		}
		return new QuotientAndRemainderResult<>(quotient, residue);
	}

	@Override
	public BigInteger euclidMeasure(ModNFE t) {
		if (t.equals(zero())) {
			return null;
		}
		return BigInteger.ZERO;
	}

	@Override
	public ModNFE projectToUnit(ModNFE t) {
//		if (primePower) {
//			return reduceComplete(complete.ringOfIntegers().projectToUnit(liftComplete(t)));
//		}
		NFE result = order.zero();
//		for (int i = 0; i < primePowerDecomposition().size(); i++) {
//			ModuloNumberFieldOrderIdeal mod = primePowerDecomposition().get(i);
//			result = order.add(mod.lift(mod.projectToUnit(mod.reduce(lift(t)))), result);
//		}
		return reduce(result);
	}

	@Override
	public Iterable<ModNFE> getUnits() {
//		getUnitGenerators();
//		return new Iterable<>() {
//			@Override
//			public Iterator<ModNFE> iterator() {
//				return new Iterator<>() {
//					private Iterator<Mod<Vector<IntE>>> it = genericModule.iterator();
//
//					@Override
//					public boolean hasNext() {
//						return it.hasNext();
//					}
//
//					@Override
//					public ModNFE next() {
//						Vector<IntE> next = genericModule.lift(it.next());
//						ModNFE result = one();
//						for (int i = 0; i < unitGenerators.size(); i++) {
//							result = multiply(power(unitGenerators.get(i), next.get(i + 1)), result);
//						}
//						return result;
//					}
//				};
//			}
//		};
		return null;
	}

	@Override
	public int krullDimension() {
		return 0;
	}

	@Override
	public IdealResult<ModNFE, ?> getIdealWithTransforms(List<ModNFE> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<ModNFE> intersect(Ideal<ModNFE> t1, Ideal<ModNFE> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<ModNFE> radical(Ideal<ModNFE> t) {
		// TODO Auto-generated method stub
		return null;
	}
}
