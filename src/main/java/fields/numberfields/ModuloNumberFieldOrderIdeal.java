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
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
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

	private NumberFieldIntegers order;
	private NumberFieldIdeal ideal;
	private RationalLattice orderAsRationalLattice;
	private SubLattice<Vector<Fraction>, Fraction, Vector<Fraction>> idealAsSubLattice;
	private boolean primePower;
	private NumberFieldIdeal primeIdeal;
	private int multiplicity;
	private CompletedNumberField complete;
	private MathMap<NFE, Ext> embedding;
	private MathMap<Ext, NFE> retraction;
	private GenericPIDModule<IntE, NFE> genericIntModule;
	private List<ModNFE> unitGenerators;
	private GenericPIDModule<IntE, Vector<IntE>> genericModule;
	private Map<ModNFE, Vector<IntE>> knownVectors;
	private List<ModuloNumberFieldOrderIdeal> primePowerDecomposition;
	private List<NFE> chineseRemainderMultipliers;

	ModuloNumberFieldOrderIdeal(NumberFieldIntegers order, NumberFieldIdeal ideal) {
		Rationals q = Rationals.q();
		this.order = order;
		this.ideal = ideal;
		List<NFE> reducedIdealBasis = new SubLattice<>(order, ideal.asSubModule(), 1.0).getModuleGenerators();
		this.orderAsRationalLattice = RationalLattice.fromIntegerLattice(order.numberField().degree(),
				new FreeModule<>(Integers.z(), order.numberField().degree()).getBasis(), 1.0);
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
		this.knownVectors = new TreeMap<>();
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(ideal);
		if (factorization.isPrimePower()) {
			this.primePower = true;
			this.primeIdeal = (NumberFieldIdeal) factorization.firstPrimeFactor();
			this.multiplicity = factorization.multiplicity(primeIdeal);
			OtherVersion<NFE, Ext, FFE, CompletedNumberField> otherVersion = order.localizeAndQuotient(primeIdeal)
					.complete(multiplicity + 1);
			this.complete = otherVersion.getField();
			this.embedding = otherVersion.getRetraction();
			this.retraction = otherVersion.getEmbedding();
		} else {
			this.primePower = false;

		}
	}

	private List<NFE> chineseRemainderMultipliers() {
		if (chineseRemainderMultipliers == null) {
			PrimaryDecompositionResult<NFE, NumberFieldIdeal> primaryDecomposition = order.primaryDecomposition(ideal);
			ChineseRemainderPreparation<NFE> crp = order
					.prepareChineseRemainderTheorem(primaryDecomposition.getPrimaryIdeals());
			this.primePowerDecomposition = new ArrayList<>();
			this.chineseRemainderMultipliers = new ArrayList<>();
			for (int i = 0; i < primaryDecomposition.getPrimaryIdeals().size(); i++) {
				primePowerDecomposition.add(primaryDecomposition.getPrimaryIdeals().get(i).modOut());
				chineseRemainderMultipliers.add(crp.getMultipliers().get(i));
			}
		}
		return chineseRemainderMultipliers;
	}

	private List<ModuloNumberFieldOrderIdeal> primePowerDecomposition() {
		chineseRemainderMultipliers();
		return primePowerDecomposition;
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
		return ideal.norm().getValue();
	}

	@Override
	public BigInteger getNumberOfUnits() throws InfinityException {
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(ideal);
		BigInteger result = BigInteger.ONE;
		for (Ideal<NFE> primeFactor : factorization.primeFactors()) {
			NumberFieldIdeal prime = (NumberFieldIdeal) primeFactor;
			BigInteger q = prime.norm().getValue();
			result = result.multiply(q.subtract(BigInteger.ONE)).multiply(q.pow(factorization.multiplicity(prime) - 1));
		}
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
		return ideal.intGenerator().getValue();
	}

	public ModNFE reduceComplete(Ext t) {
		return reduce(retraction.evaluate(complete.round(t, multiplicity)));// new ModNFE(complete.round(t,
																			// multiplicity));
	}

	public ModNFE reduce(NFE t) {
		Rationals q = Rationals.q();
		NFE reduced = genericIntModule.lift(genericIntModule.reduce(t));
		Vector<Fraction> asFractions = Vector.mapVector(q.getEmbeddingMap(), order.asVector(reduced));
		Vector<Fraction>	closestLatticePoint = orderAsRationalLattice.getVectorSpace().closestLatticePoint(asFractions, idealAsSubLattice);
		Vector<Fraction> diff = orderAsRationalLattice.getVectorSpace().subtract(asFractions, closestLatticePoint);
		return new ModNFE(order.fromVector(Vector.mapVector(q.getAsIntegerMap(), diff)));
	}

	public Ext liftComplete(ModNFE t) {
		return complete.round(embedding.evaluate(lift(t)), multiplicity);
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
		// return reduceComplete(complete.multiply(liftComplete(t1), liftComplete(t2)));
	}

	@Override
	public boolean isUnit(ModNFE t) {
		if (primePower) {
			return complete.valuation(liftComplete(t)).equals(Value.ZERO);
		}
		for (ModuloNumberFieldOrderIdeal mod : primePowerDecomposition()) {
			if (!mod.isUnit(mod.reduce(lift(t)))) {
				return false;
			}
		}
		return true;
	}

	@Override
	public ModNFE inverse(ModNFE t) {
		if (primePower) {
			return reduceComplete(complete.inverse(liftComplete(t)));
		}
		NFE result = order.zero();
		for (int i = 0; i < primePowerDecomposition().size(); i++) {
			ModuloNumberFieldOrderIdeal mod = primePowerDecomposition().get(i);
			result = order.add(mod.lift(mod.inverse(mod.reduce(lift(t)))), result);
		}
		return reduce(result);
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

	@Override
	public boolean isDivisible(ModNFE dividend, ModNFE divisor) {
		if (divisor.equals(zero())) {
			return false;
		}
		if (primePower) {
			return complete.valuation(liftComplete(dividend)).compareTo(complete.valuation(liftComplete(divisor))) >= 0;
		}
		for (ModuloNumberFieldOrderIdeal mod : primePowerDecomposition()) {
			if (!mod.isDivisible(mod.reduce(lift(dividend)), mod.reduce(lift(divisor)))) {
				return false;
			}
		}
		return true;
		// return order.add(ideal,
		// order.getIdeal(Collections.singletonList(divisor.lift))).contains(dividend.lift);
	}

	@Override
	public QuotientAndRemainderResult<ModNFE> quotientAndRemainder(ModNFE dividend, ModNFE divisor) {
		if (primePower) {
			QuotientAndRemainderResult<Ext> qr = complete.ringOfIntegers().quotientAndRemainder(liftComplete(dividend),
					liftComplete(divisor));
			return new QuotientAndRemainderResult<>(reduceComplete(qr.getQuotient()),
					reduceComplete(qr.getRemainder()));
		}
		NFE quotient = order.zero();
		NFE remainder = order.zero();
		for (int i = 0; i < primePowerDecomposition().size(); i++) {
			ModuloNumberFieldOrderIdeal mod = primePowerDecomposition().get(i);
			QuotientAndRemainderResult<ModNFE> qr = mod.quotientAndRemainder(mod.reduce(lift(dividend)),
					mod.reduce(lift(divisor)));
			quotient = order.add(mod.lift(qr.getQuotient()), quotient);
			remainder = order.add(mod.lift(qr.getRemainder()), remainder);
		}
		return new QuotientAndRemainderResult<>(reduce(quotient), reduce(remainder));
		// List<NFE> generators = new ArrayList<>();
//		generators.add(divisor.lift);
//		generators.addAll(ideal.generators());
//		IdealResult<NFE, NumberFieldIdeal> asIdeal = order.getIdealWithTransforms(generators);
//		List<NFE> generate = asIdeal.getIdeal().generate(dividend.lift);
//		ModNFE residue = reduce(asIdeal.getIdeal().residue(dividend.lift));
//		ModNFE quotient = zero();
//		for (int i = 0; i < generate.size(); i++) {
//			quotient = add(quotient,
//					reduce(order.multiply(generate.get(i), asIdeal.getGeneratorExpressions().get(i).get(0))));
//		}
//		return new QuotientAndRemainderResult<>(quotient, residue);
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
		if (primePower) {
			return reduceComplete(complete.ringOfIntegers().projectToUnit(liftComplete(t)));
		}
		NFE result = order.zero();
		for (int i = 0; i < primePowerDecomposition().size(); i++) {
			ModuloNumberFieldOrderIdeal mod = primePowerDecomposition().get(i);
			result = order.add(mod.lift(mod.projectToUnit(mod.reduce(lift(t)))), result);
		}
		return reduce(result);
	}

	@Override
	public Iterable<ModNFE> getUnits() {
		getUnitGenerators();
		return new Iterable<>() {
			@Override
			public Iterator<ModNFE> iterator() {
				return new Iterator<>() {
					private Iterator<Mod<Vector<IntE>>> it = genericModule.iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public ModNFE next() {
						Vector<IntE> next = genericModule.lift(it.next());
						ModNFE result = one();
						for (int i = 0; i < unitGenerators.size(); i++) {
							result = multiply(power(unitGenerators.get(i), next.get(i + 1)), result);
						}
						return result;
					}
				};
			}
		};
	}

	public Vector<IntE> asUnitGeneratorVector(ModNFE t) {
		if (!isUnit(t)) {
			throw new ArithmeticException(t + " is not a unit!");
		}
		while (!knownVectors.containsKey(t)) {
			FreeModule<IntE> free = new FreeModule<>(Integers.z(), getUnitGenerators().size());
			Map<ModNFE, Vector<IntE>> giantSteps = new TreeMap<>();
			Vector<IntE> asVector = genericModule.lift(genericModule.getRandomElement());
			ModNFE babyStep = fromUnitGeneratorVector(asVector);
			if (!knownVectors.containsKey(babyStep)) {
				knownVectors.put(babyStep, asVector);
			}
			ModNFE giantStep = divideChecked(t, babyStep);
			if (!giantSteps.containsKey(giantStep)) {
				giantSteps.put(giantStep, asVector);
			}
			if (knownVectors.containsKey(giantStep)) {
				Vector<IntE> result = genericModule
						.lift(genericModule.reduce(free.add(knownVectors.get(giantStep), asVector)));
				knownVectors.put(t, result);
			} else if (giantSteps.containsKey(babyStep)) {
				Vector<IntE> result = genericModule
						.lift(genericModule.reduce(free.add(giantSteps.get(babyStep), asVector)));
				knownVectors.put(t, result);
			}
		}
		return knownVectors.get(t);
	}

	public ModNFE fromUnitGeneratorVector(Vector<IntE> t) {
		ModNFE result = one();
		List<ModNFE> unitGenerators = getUnitGenerators();
		for (int i = 0; i < unitGenerators.size(); i++) {
			result = multiply(power(unitGenerators.get(i), t.get(i + 1)), result);
		}
		return result;
	}

	public List<ModNFE> getUnitGenerators() {
		if (primePower) {
			return getUnitGeneratorsPrimePower();
		}
		if (unitGenerators == null) {
			Integers z = Integers.z();
			List<ModNFE> modUnitGenerators = new ArrayList<>();
			List<List<IntE>> syzygies = new ArrayList<>();
			int numGenerators = 0;
			ModNFE multiplierSum = zero();
			for (int i = 0; i < primePowerDecomposition().size(); i++) {
				List<ModNFE> generators = primePowerDecomposition().get(i).getUnitGenerators();
				ModNFE multiplier = reduce(chineseRemainderMultipliers().get(i));
				for (int j = 0; j < modUnitGenerators.size(); j++) {
					modUnitGenerators.set(j, add(multiplier, modUnitGenerators.get(j)));
				}
				for (ModNFE generator : generators) {
					modUnitGenerators.add(add(multiplierSum, multiply(multiplier, generator)));
				}
				for (List<IntE> syzygy : syzygies) {
					for (int j = 0; j < generators.size(); j++) {
						syzygy.add(z.zero());
					}
				}
				for (List<IntE> syzygy : primePowerDecomposition().get(i).genericModule.getModuleGeneratorRelations()) {
					List<IntE> augmented = new ArrayList<>();
					for (int j = 0; j < numGenerators; j++) {
						augmented.add(z.zero());
					}
					augmented.addAll(syzygy);
					syzygies.add(augmented);
				}
				numGenerators += generators.size();
				multiplierSum = add(multiplier, multiplierSum);
			}
			FreeModule<IntE> free = new FreeModule<>(z, numGenerators);
			GenericPIDModule<IntE, Vector<IntE>> asGenericModule = GenericPIDModule.fromSyzygies(free, syzygies);
			unitGenerators = new ArrayList<>();
			List<List<IntE>> unitSyzygies = new ArrayList<>();
			int nonTrivialIndex = -1;
			for (int i = 0; i < asGenericModule.diagonalBasis().size(); i++) {
				if (z.isUnit(asGenericModule.diagonalRanks().get(i))) {
					continue;
				}
				if (nonTrivialIndex < 0) {
					nonTrivialIndex = i;
					free = new FreeModule<>(z, asGenericModule.diagonalBasis().size() - nonTrivialIndex);
				}
				Vector<IntE> basis = asGenericModule.lift(asGenericModule.diagonalBasis().get(i));
				ModNFE asElement = one();
				for (int j = 0; j < modUnitGenerators.size(); j++) {
					asElement = multiply(power(modUnitGenerators.get(j), basis.get(j + 1)), asElement);
				}
				unitGenerators.add(asElement);
				unitSyzygies.add(free.scalarMultiply(asGenericModule.diagonalRanks().get(i),
						free.getUnitVector(i - nonTrivialIndex + 1)).asList());
			}
			genericModule = GenericPIDModule.fromSyzygies(free, unitSyzygies);

		}
		return unitGenerators;
	}

	private List<ModNFE> getUnitGeneratorsPrimePower() {
		if (unitGenerators == null) {
			Integers z = Integers.z();
			FiniteField field = complete.ringOfIntegers().reduction();
			FFE primitiveRoot = field.primitiveRoot();
			UnivariatePolynomialRing<Ext> polynomialRing = complete.getUnivariatePolynomialRing();
			UnivariatePolynomial<Ext> unity = polynomialRing.toUnivariate(polynomialRing.subtract(
					polynomialRing.getVarPower(field.getNumberOfElements().intValueExact() - 1), polynomialRing.one()));
			Ext lifted = complete.ringOfIntegers().henselLift(unity, primitiveRoot, multiplicity);
			List<ModNFE> generators = new ArrayList<>();
			generators.add(reduceComplete(lifted));
			for (int k = 1; k < multiplicity; k++) {
				NumberFieldIdeal idealPower = order.power(primeIdeal, k);
				for (NFE generator : idealPower.asSubModule().getModuleGenerators()) {
					generators.add(reduce(order.add(generator, order.one())));
				}
			}
			Group<ModNFE> multiplicativeGroup = getMultiplicativeGroup();
			FreeModule<IntE> free = new FreeModule<>(z, generators.size());
			List<Vector<IntE>> syzygies = new ArrayList<>();
			Map<ModNFE, Vector<IntE>> knownVectors = new TreeMap<>();
			knownVectors.put(one(), free.zero());
			syzygies.add(free.scalarMultiply(field.getNumberOfUnits(), free.getUnitVector(1)));
			knownVectors.put(generators.get(0), free.getUnitVector(1));
			for (int i = 1; i < generators.size(); i++) {
				if (knownVectors.containsKey(generators.get(i))) {
					syzygies.add(free.subtract(knownVectors.get(generators.get(i)), free.getUnitVector(i + 1)));
				} else {
					IntE generatorOrder = z.getInteger(multiplicativeGroup.getOrder(generators.get(i)));
					knownVectors.put(generators.get(i), free.getUnitVector(i + 1));
					syzygies.add(free.scalarMultiply(generatorOrder, free.getUnitVector(i + 1)));
				}
			}
			FreeSubModule<IntE, Vector<IntE>> syzygyModule = new FreeSubModule<>(free, syzygies);
			GenericPIDModule<IntE, Vector<IntE>> asGenericModule = new GenericPIDModule<>(free, syzygyModule);
			BigInteger numberOfUnits = getNumberOfUnits();
			while (!numberOfUnits.equals(asGenericModule.getNumberOfElements())) {
				Vector<IntE> random = asGenericModule.lift(asGenericModule.getRandomElement());
				ModNFE element = one();
				for (int i = 0; i < generators.size(); i++) {
					IntE power = random.get(i + 1);
					element = multiply(power(generators.get(i), power), element);
				}
				if (!knownVectors.containsKey(element)) {
					knownVectors.put(element, random);
				} else {
					Vector<IntE> newSyzygy = free.subtract(knownVectors.get(element), random);
					if (!syzygyModule.contains(newSyzygy)) {
						syzygies.add(newSyzygy);
						syzygyModule = new FreeSubModule<>(free, syzygies);
						asGenericModule = new GenericPIDModule<>(free, syzygyModule);
					}
				}
			}
			unitGenerators = new ArrayList<>();
			List<List<IntE>> unitSyzygies = new ArrayList<>();
			int nonTrivialIndex = -1;
			for (int i = 0; i < asGenericModule.diagonalBasis().size(); i++) {
				if (z.isUnit(asGenericModule.diagonalRanks().get(i))) {
					continue;
				}
				if (nonTrivialIndex < 0) {
					nonTrivialIndex = i;
					free = new FreeModule<>(z, asGenericModule.diagonalBasis().size() - nonTrivialIndex);
				}
				Vector<IntE> basis = asGenericModule.lift(asGenericModule.diagonalBasis().get(i));
				ModNFE asElement = one();
				for (int j = 0; j < generators.size(); j++) {
					asElement = multiply(power(generators.get(j), basis.get(j + 1)), asElement);
				}
				unitGenerators.add(asElement);
				unitSyzygies.add(free.scalarMultiply(asGenericModule.diagonalRanks().get(i),
						free.getUnitVector(i - nonTrivialIndex + 1)).asList());
			}
			genericModule = GenericPIDModule.fromSyzygies(free, unitSyzygies);
		}
		return unitGenerators;
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
