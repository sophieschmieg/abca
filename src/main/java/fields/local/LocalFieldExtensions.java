package fields.local;
//package fields.local;
//
//import java.math.BigInteger;
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.List;
//import java.util.Map;
//import java.util.TreeMap;
//
//import fields.helper.AbstractElement;
//import fields.helper.AbstractFieldExtension;
//import fields.helper.FieldAutomorphism;
//import fields.helper.FieldEmbedding;
//import fields.helper.GaloisGroup;
//import fields.helper.GenericExtensionField;
//import fields.helper.GenericExtensionField.GenericExtensionFieldElement;
//import fields.interfaces.Element;
//import fields.interfaces.ExtensionFieldElement;
//import fields.interfaces.Field;
//import fields.interfaces.FieldExtension;
//import fields.interfaces.LocalField;
//import fields.interfaces.LocalRing;
//import fields.interfaces.MathMap;
//import fields.interfaces.Polynomial;
//import fields.local.LocalFieldExtension.Ext;
//import fields.polynomials.MultivariatePolynomialRing;
//import fields.polynomials.UnivariatePolynomial;
//import fields.polynomials.UnivariatePolynomialRing;
//import fields.vectors.FreeModule;
//import fields.vectors.Matrix;
//import fields.vectors.MatrixAlgebra;
//import fields.vectors.Vector;
//
//public class LocalFieldExtension<Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>>
//		extends
//		AbstractFieldExtension<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>>
//		implements LocalField<Ext<Base>, ReducedExt>,
//		FieldExtension<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> {
//	//private LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base;
//	private LocalField<Base, Reduced> prime;
//	private ReductionExtensionType reduction;
//	private int ramificationIndex;
//	private int residueDegree;
////	private int fullRamificationIndex;
////	private int fullResidueDegree;
//	private Ext<Base> uniformizer;
//	private Ext<Base> liftOfReductionGenerator;
//	private List<Ext<Base>> reductionGeneratorUniformizerBase;
//	private Matrix<Base> gammaToReductionUniformizer;
////	private boolean trivial;
//	private LocalRing<Ext<Base>, ReducedExt> localRing;
//	private Map<Ext<Base>, List<Ext<Base>>> conjugates;
//	private Map<UnivariatePolynomial<Ext<Base>>, FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>>> extensions;
//	private LocalField<Base, Reduced> highAccuracyPrime;
//	private MatrixAlgebra<Base> highAccuracyAlgebra;
//	private UnivariatePolynomialRing<Base> highAccuracyPolynomials;
////	private LocalField<Base, Reduced> base;
////	private UnramifiedExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramified;
////	private TotallyRamifiedExtension<Base, Reduced> totallyRamified;
////	private MultivariatePolynomialRing<Base> ring;
////	private Ideal<Polynomial<Base>> ideal;
////	private Polynomial<Base> unramifiedMinimalPolynomial;
////	private Polynomial<Base> totallyRamifiedMinimalPolynomial;
////	private Map<Monomial, Polynomial<Base>> powers;
////	private MatrixAlgebra<Base> baseAlgebra;
////	private Matrix<Base> baseChange;
////	private Ext<Base> zero;
////	private Ext<Base> one;
////	private Ext<Base> reductionGenerator;
////	private Ext<Base> uniformizer;
//
//	public static class Ext<Base extends Element<Base>> extends AbstractElement<Ext<Base>>
//			implements ExtensionFieldElement<Base, Ext<Base>> {
//		private UnivariatePolynomial<Base> e;
//		private Value valuation;
////		private Field<Base> base;
////		private int ramificationIndex;
////		private int residueDegree;
////		private MatrixAlgebra<Base> baseAlgebra;
////		private Matrix<Base> baseChange;
//
//		private Ext(UnivariatePolynomial<Base> e) {
//			this.e = e;
//		}
////		private Ext(Field<Base> base, int residueDegree, int ramificationIndex, MatrixAlgebra<Base> baseAlgebra,
////				Matrix<Base> baseChange, Polynomial<Base> e) {
////			this.e = e;
////			this.base = base;
////			this.residueDegree = residueDegree;
////			this.ramificationIndex = ramificationIndex;
////			this.baseAlgebra = baseAlgebra;
////			this.baseChange = baseChange;
////		}
//
//		@Override
//		public int compareTo(Ext<Base> o) {
//			return e.compareTo(o.e);
//		}
//		//
////		@Override
////		public String toString() {
////			return e.toString(new String[] { "α", "τ" }, true);
////		}
//
//		@Override
//		public String toString() {
//			return e.toString("α", true);
//		}
//
//		@Override
//		public UnivariatePolynomial<Base> asPolynomial() {
//			return e;
////			Vector<Base> asVector = ((MultivariatePolynomialRing<Base>) e.getPolynomialRing()).asVector(e,
////					new int[] { residueDegree, ramificationIndex });
////			Vector<Base> asBaseChangeVector = baseAlgebra.multiply(baseChange, asVector);
////			return base.getUnivariatePolynomialRing().getPolynomial(asBaseChangeVector.asList());
//		}
//
//	}
//
//	LocalFieldExtension(LocalField<Base, Reduced> base, ReductionExtensionType reductionAsTrivialExtension) {
//		super(base);
//		//this.base = this;
//		this.reduction = reductionAsTrivialExtension
//				.makeExtension(base.reduction().getUnivariatePolynomialRing().getVar());
//		this.ramificationIndex = 1;
//		this.residueDegree = 1;
//		this.uniformizer = getEmbedding(base.uniformizer());
//		this.liftOfReductionGenerator = zero();
//		this.reductionGeneratorUniformizerBase = Collections.singletonList(one());
//		this.gammaToReductionUniformizer = matrixAlgebra().one();
//		this.localRing = new LocalRingImplementation<>(this, base.ringOfIntegers().toString());
//		// this.trivial = true;
//		this.prime = base;
//		this.conjugates = new TreeMap<>();
//		this.extensions = new TreeMap<>();
//		this.highAccuracyPrime = prime;
//		this.highAccuracyPolynomials = prime.getUnivariatePolynomialRing();
//		this.highAccuracyAlgebra = matrixAlgebra();
//	}
////	private static <Base extends Element<Base>, Ext extends ExtensionFieldElement<Base, Ext>, T extends ExtensionFieldElement<Base, T>> GenericExtensionField<?> asGenericExtensionField(GenericExtensionField<Base> extension, UnivariatePolynomial<T> minimalPolynomial) {
////	GenericExtensionField<Base> field = new GenericExtensionField<>(extension.getBaseField());
////	//.getExtension(minimalPolynomial).extension();
////	}
////	
////	private static <Base extends Element<Base>, Ext extends ExtensionFieldElement<Base, Ext>> GenericExtensionField<?> asGenericExtensionField(FieldExtension<Base, Ext> extension) {
////		if (extension.getBaseField() instanceof FieldExtension<?, ?>) {
////			GenericExtensionField<?> field = asGenericExtensionField((FieldExtension<?, ?>)extension.getBaseField());
////			UnivariatePolynomial<?> genericMinimalPolynomial = field.getUnivariatePolynomialRing().getEmbedding( extension.minimalPolynomial(), new MathMap<>() {
////				@Override
////				public GenericExtensionFieldElement<?> evaluate(Base t) {
////					return ExtensionFieldElement<?>;
////				}});
////		}
////		return new GenericExtensionField<>(extension.minimalPolynomial(), extension.getBaseField());
////	}
////	
////	@SuppressWarnings("unchecked")
////	static <Base extends Element<Base>, Ext extends ExtensionFieldElement<Base, Ext>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>> LocalFieldExtension<?, ?, ?, ?> createLocalFieldExtension(
////			UnivariatePolynomial<Ext> minimalPolynomial, LocalFieldExtension<Base, Ext, Reduced, ReducedExt> base) {
////		GenericExtensionField<?> baseAsGeneric;
////		LocalFieldExtension<?, ?, ?, ?> currentBase = base;
////		FieldExtension<?, ?> extAsGeneric =
////		baseAsGeneric.getExtension(baseAsGeneric.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial, new MathMap<>() {
////					@Override
////					public GenericExtensionFieldElement<Base> evaluate( t) {
////						return baseAsGeneric.fromPolynomial(t.asPolynomial());
////					}})).extension();
////
////		while (currentBase.getBaseField() instanceof LocalFieldExtension<?, ?, ?, ?>) {
////			currentBase = (LocalFieldExtension<?, ?, ?, ?>)currentBase.getBaseField();
////		}
////		 = new GenericExtensionField<>(base.minimalPolynomial(), base.getBaseField());
////				return createLocalFieldExtension(extAsGeneric.minimalPolynomial(), base.getBaseField(), (FieldExtension<Reduced, ReducedExt>)base.reduction().getExtension(base.reduction().getUnivariatePolynomialRing().getVar()));
////	}
////	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> Extension<Base, Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> createLocalFieldExtension(
////			UnivariatePolynomial<Base> minimalPolynomial, LocalField<Base, Reduced> field,
////			ReductionExtensionType reductionAsTrivialExtension) {
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base = new LocalFieldExtension<>(field,
////				reductionAsTrivialExtension);
////		UnivariatePolynomial<Ext<Base>> embeddedMinimalPolynomial = base.getUnivariatePolynomialRing()
////				.getEmbedding(minimalPolynomial, base.getEmbeddingMap());
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> extension = createLocalFieldExtension(
////				embeddedMinimalPolynomial, base).getField();
////		return new Extension<>(extension, extension.getEmbeddingMap(), extension.asBaseFieldElementMap());
////	}
////
////	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> Extension<URExt<Base>, Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> createLocalFieldExtension(
////			UnivariatePolynomial<URExt<Base>> minimalPolynomial,
////			UnramifiedExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramified) {
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base = new LocalFieldExtension<>(
////				unramified);
////		UnivariatePolynomial<Ext<Base>> embeddedMinimalPolynomial = base.getUnivariatePolynomialRing()
////				.getEmbedding(minimalPolynomial, new MathMap<>() {
////					@Override
////					public Ext<Base> evaluate(URExt<Base> t) {
////						return t.asExt(base);
////					}
////				});
////		FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> extension = createLocalFieldExtension(
////				embeddedMinimalPolynomial, base);
////		return new Extension<URExt<Base>, Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>>(
////				extension.getField(), new MathMap<>() {
////					@Override
////					public Ext<Base> evaluate(URExt<Base> t) {
////						return extension.getEmbedding(t.asExt(base));
////					}
////				}, new MathMap<>() {
////					@Override
////					public URExt<Base> evaluate(Ext<Base> t) {
////						return unramified.fromPolynomial(extension.asVector(t).get(1).asPolynomial());
////					}
////				});
////	}
////
////	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> Extension<TRExt<Base>, Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> createLocalFieldExtension(
////			UnivariatePolynomial<TRExt<Base>> minimalPolynomial,
////			TotallyRamifiedExtension<Base, Reduced> totallyRamified,
////			ReductionExtensionType reductionAsTrivialExtension) {
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base = new LocalFieldExtension<>(
////				totallyRamified, reductionAsTrivialExtension);
////		UnivariatePolynomial<Ext<Base>> embeddedMinimalPolynomial = base.getUnivariatePolynomialRing()
////				.getEmbedding(minimalPolynomial, new MathMap<>() {
////					@Override
////					public Ext<Base> evaluate(TRExt<Base> t) {
////						return t.asExt(base);
////					}
////				});
////		FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> extension = createLocalFieldExtension(
////				embeddedMinimalPolynomial, base);
////		return new Extension<TRExt<Base>, Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>>(
////				extension.getField(), new MathMap<>() {
////					@Override
////					public Ext<Base> evaluate(TRExt<Base> t) {
////						return extension.getEmbedding(t.asExt(base));
////					}
////				}, new MathMap<>() {
////					@Override
////					public TRExt<Base> evaluate(Ext<Base> t) {
////						return totallyRamified.fromPolynomial(extension.asVector(t).get(1).asPolynomial());
////					}
////				});
////	}
////	
////	private static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> createUnramifiedExtension(
////			UnivariatePolynomial<Ext<Base>> minimalPolynomial, UnivariatePolynomial<ReducedExt> reductionMinimalPolynomial,
////			LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base) {
////		FieldEmbedding<Base, URExt<Base>, UnramifiedExtension<Base, Reduced, ReducedExt,ReductionExtensionType>> unramifiedExtension = base.unramified.getEmbeddedExtension(base.unramified.ringOfIntegers().liftUnivariatePolynomial(reductionMinimalPolynomial));
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> extension = new LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>(
////				base.getBaseField(), unramifiedExtension.getField(), base.totallyRamified);
////	return new FieldEmbedding<>(extension, base, new MathMap<>() {
////		@Override
////		public Ext<Base> evaluate(Ext<Base> t) {
////			return extension.fromMultivariatePolynomial();
////		}});
////	}
////	
////	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> createLocalFieldExtension(
////			UnivariatePolynomial<Ext<Base>> minimalPolynomial,
////			LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base) {
////		if (minimalPolynomial.degree() == 1) {
////			return new FieldEmbedding<>(base);
////		}
////		if (!base.isIrreducible(minimalPolynomial)) {
////			throw new ArithmeticException("Polynomial not irreducible!");
////		}
////		LocalRing<Ext<Base>, ReducedExt> baseRing = base.ringOfIntegers();
////		UnivariatePolynomialRing<ReducedExt> baseReducedPolynomials = base.reduction().getUnivariatePolynomialRing();
////		minimalPolynomial = baseRing.integralMinimalPolynomial(minimalPolynomial);
////		UnivariatePolynomial<ReducedExt> baseReduced = baseReducedPolynomials.toUnivariate(
////				baseReducedPolynomials.squareFreeFactorization(baseRing.reduceUnivariatePolynomial(minimalPolynomial))
////						.keySet().iterator().next());
////		if (baseReduced.degree() == minimalPolynomial.degree()) {
////		return createUnramifiedExtension(minimalPolynomial, baseReduced, base);
////		}
////		if (baseRing.isEisenstein(minimalPolynomial)) {
////			return new LocalFieldExtension<>(new TotallyRamifiedExtension<>(minimalPolynomial, base), base.reduction());
////		}
////		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramifiedExtension = base;
////		LocalRing<Ext<Base>, ReducedExt> unramifiedExtensionIntegers = unramifiedExtension.ringOfIntegers();
////		UnivariatePolynomialRing<Ext<Base>> unramifiedPolynomials = unramifiedExtension.getUnivariatePolynomialRing();
////		UnivariatePolynomial<Ext<Base>> minimalPolynomialOverUnramifiedExtension = minimalPolynomial;
////		UnivariatePolynomialRing<ReducedExt> reducedPolynomials = unramifiedExtension.reduction()
////				.getUnivariatePolynomialRing();
////		UnivariatePolynomial<ReducedExt> reduced = baseReduced;
////		while (true) {
////			if (minimalPolynomialOverUnramifiedExtension.degree() == 1) {
////				return unramifiedExtension;
////			} else if (reduced.degree() != 1) {
////				FieldEmbedding<Base, LocalExtensionFieldElement<Base>> extension = base.unramified.getEmbeddedExtension(
////						unramifiedExtension.unramified.ringOfIntegers().liftUnivariatePolynomial(reduced));
////				unramifiedExtension = new LocalFieldExtension<>(base.getBaseField(),
////						(UnramifiedExtension<Base, Reduced, ReducedExt>) extension.getField(), base.totallyRamified);
////				unramifiedExtensionIntegers = unramifiedExtension.ringOfIntegers();
////				unramifiedPolynomials = unramifiedExtension.getUnivariatePolynomialRing();
////				minimalPolynomialOverUnramifiedExtension = unramifiedPolynomials
////						.toUnivariate(unramifiedExtension
////								.factorization(unramifiedPolynomials.getEmbedding(minimalPolynomial,
////										unramifiedExtension.getEmbeddingMap()))
////								.getFactors().keySet().iterator().next());
////				reducedPolynomials = unramifiedExtension.reduction().getUnivariatePolynomialRing();
////				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
////						.squareFreeFactorization(unramifiedExtensionIntegers
////								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
////						.keySet().iterator().next());
////			} else if (!reduced.univariateCoefficient(0).equals(unramifiedExtension.reduction().zero())) {
////				LocalExtensionFieldElement<Base> constant = unramifiedExtension.lift(reduced.univariateCoefficient(0));
////				Polynomial<LocalExtensionFieldElement<Base>> transform = unramifiedPolynomials
////						.getPolynomial(unramifiedExtension.negative(constant), unramifiedExtension.one());
////				minimalPolynomialOverUnramifiedExtension = unramifiedPolynomials
////						.substitute(minimalPolynomialOverUnramifiedExtension, Collections.singletonList(transform));
////				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
////						.squareFreeFactorization(unramifiedExtensionIntegers
////								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
////						.keySet().iterator().next());
////			} else if (unramifiedExtension.valuation(minimalPolynomialOverUnramifiedExtension.univariateCoefficient(0))
////					.equals(new Value(1))) {
////				return new LocalFieldExtension<>(unramifiedExtension,
////						new TotallyRamifiedExtension<>(minimalPolynomialOverUnramifiedExtension, unramifiedExtension));
////			} else {
////				minimalPolynomialOverUnramifiedExtension = unramifiedExtensionIntegers
////						.integralMinimalPolynomial(minimalPolynomialOverUnramifiedExtension);
////				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
////						.squareFreeFactorization(unramifiedExtensionIntegers
////								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
////						.keySet().iterator().next());
////				if (unramifiedExtension.valuation(minimalPolynomialOverUnramifiedExtension.univariateCoefficient(0))
////						.compareTo(new Value(1)) > 0) {
////					GenericExtensionField<LocalExtensionFieldElement<Base>> extension = new GenericExtensionField<>(
////							minimalPolynomialOverUnramifiedExtension, unramifiedExtension);
////					int ramificationIndex = minimalPolynomialOverUnramifiedExtension.degree() / reduced.degree();
////					GenericExtensionFieldElement<LocalExtensionFieldElement<Base>> alpha = extension.alpha();
////					int alphaValuation = unramifiedExtension.valuation(extension.norm(alpha)).value();
////					while (true) {
////						if (BigInteger.valueOf(alphaValuation).gcd(BigInteger.valueOf(ramificationIndex))
////								.intValueExact() == 1) {
////							int power = BigInteger.valueOf(alphaValuation)
////									.modInverse(BigInteger.valueOf(ramificationIndex)).intValueExact();
////							int ramificationPower = (power * alphaValuation - 1) / ramificationIndex;
////							minimalPolynomialOverUnramifiedExtension = extension
////									.minimalPolynomial(extension.divide(extension.power(alpha, power),
////											extension.power(extension.getEmbedding(unramifiedExtension.uniformizer()),
////													ramificationPower)));// ,
////							reduced = reducedPolynomials.toUnivariate(reducedPolynomials
////									.squareFreeFactorization(unramifiedExtensionIntegers
////											.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
////									.keySet().iterator().next());
////							break;
////						} else {
////							int gcd = BigInteger.valueOf(alphaValuation).gcd(BigInteger.valueOf(ramificationIndex))
////									.intValueExact();
////							alpha = extension.add(extension.power(alpha, ramificationIndex / gcd),
////									extension.getEmbedding(unramifiedExtension.power(unramifiedExtension.uniformizer(),
////											alphaValuation / gcd)));
////							alphaValuation = unramifiedExtension.valuation(extension.norm(alpha)).value();
////						}
////					}
////				}
////			}
////		}
////	}
////
////	private LocalFieldExtension(LocalField<Base, Reduced> base, ReductionExtensionType reductionAsTrivialExtension) {
////		super(base.getUnivariatePolynomialRing().getVar(), base);
////		this.unramified = new UnramifiedExtension<>(base.getUnivariatePolynomialRing().getVar(), base,
////				reductionAsTrivialExtension);
////		this.totallyRamified = new TotallyRamifiedExtension<>(base.getUnivariatePolynomialRing().getVar(), base);
////		init();
////	}
////
////	private LocalFieldExtension(UnramifiedExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramified) {
////		super(unramified.minimalPolynomial(), unramified.getBaseField());
////		this.unramified = unramified;
////		this.totallyRamified = new TotallyRamifiedExtension<>(
////				unramified.getBaseField().getUnivariatePolynomialRing().getVar(), unramified.getBaseField());
////		init();
////	}
////
////	private LocalFieldExtension(TotallyRamifiedExtension<Base, Reduced> totallyRamified,
////			ReductionExtensionType reductionAsTrivialExtension) {
////		super(totallyRamified.minimalPolynomial(), totallyRamified.getBaseField());
////		LocalField<Base, Reduced> base = totallyRamified.getBaseField();
////		this.unramified = new UnramifiedExtension<>(base.getUnivariatePolynomialRing().getVar(), base,
////				reductionAsTrivialExtension);
////		this.totallyRamified = totallyRamified;
////		init();
////	}
////
////	private static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> UnivariatePolynomial<Base> combinedMinimalPolynomial(
////			LocalField<Base, Reduced> base,
////			UnramifiedExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramified,
////			TotallyRamifiedExtension<Base, Reduced> totallyRamified) {
////		MultivariatePolynomialRing<Base> ring = (MultivariatePolynomialRing<Base>) AbstractPolynomialRing
////				.getPolynomialRing(base, 2, Monomial.GREVLEX);
////		List<Polynomial<Base>> minimalPolynomials = new ArrayList<>();
////		minimalPolynomials.add(ring.getEmbedding(unramified.minimalPolynomial(), new int[] { 0 }));
////		minimalPolynomials.add(ring.getEmbedding(totallyRamified.minimalPolynomial(), new int[] { 1 }));
////		Ideal<Polynomial<Base>> ideal = ring.getIdeal(minimalPolynomials);
////		CoordinateRing<Base> cr = new CoordinateRing<Base>(ring, ideal);
////		return AbstractFieldExtension.minimalPolynomial(cr.getEmbedding(ring.add(ring.getVar(1), ring.getVar(2))),
////				unramified.degree() * totallyRamified.degree(), cr, base, new MathMap<>() {
////					@Override
////					public Vector<Base> evaluate(CoordinateRingElement<Base> t) {
////						return ring.asVector(t.getElement(),
////								new int[] { unramified.degree(), totallyRamified.degree() });
////					}
////				});
////	}
////
////	private LocalFieldExtension(LocalField<Base, Reduced> base,
////			UnramifiedExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramified,
////			TotallyRamifiedExtension<Base, Reduced> totallyRamified) {
////		super(combinedMinimalPolynomial(base, unramified, totallyRamified), base);
////		this.unramified = unramified;
////		this.totallyRamified = totallyRamified;
////		init();
////	}
////
////	private void init() {
////		this.base = unramified.getBaseField();
////		this.baseAlgebra = new FiniteVectorSpace<>(base, unramified.degree() * totallyRamified.degree())
////				.matrixAlgebra();
////		this.ring = (MultivariatePolynomialRing<Base>) AbstractPolynomialRing.getPolynomialRing(base, 2,
////				Monomial.GREVLEX);
////		this.unramifiedMinimalPolynomial = this.ring.getEmbedding(unramified.minimalPolynomial(), new int[] { 0 });
////		this.totallyRamifiedMinimalPolynomial = this.ring.getEmbedding(totallyRamified.minimalPolynomial(),
////				new int[] { 1 });
////		this.zero = fromMultivariatePolynomialNoReduce(ring.zero());
////		this.one = fromMultivariatePolynomialNoReduce(ring.one());
////		List<Polynomial<Base>> idealGen = new ArrayList<>();
////		idealGen.add(unramifiedMinimalPolynomial);
////		idealGen.add(totallyRamifiedMinimalPolynomial);
////		this.ideal = this.ring.getIdeal(idealGen);
////		this.powers = new TreeMap<>();
////		for (int i = 0; i < 2 * unramified.degree(); i++) {
////			for (int j = 0; j < 2 * totallyRamified.degree(); j++) {
////				Monomial m = ring.getMonomial(new int[] { i, j });
////				this.powers.put(m, this.ideal.residue(ring.getEmbedding(base.one(), m)));
////			}
////		}
////		this.reductionGenerator = fromMultivariatePolynomial(ring.getVar(1));
////		if (totallyRamified.degree() > 1) {
////			this.uniformizer = fromMultivariatePolynomialNoReduce(ring.getVar(2));
////		} else {
////			this.uniformizer = fromMultivariatePolynomialNoReduce(ring.getEmbedding(base.uniformizer()));
////		}
////		if (unramified.degree() == 1 || totallyRamified.degree() == 1) {
////			this.baseChange = this.baseAlgebra.one();
////		} else {
////			Ext<Base> power = one;
////			Ext<Base> alpha = alpha();
////			List<Vector<Base>> generators = new ArrayList<>();
////			for (int i = 0; i < unramified.degree() * totallyRamified.degree(); i++) {
////				generators.add(ring.asVector(power.e, new int[] { unramified.degree(), totallyRamified.degree() }));
////				power = multiply(power, alpha);
////			}
////			this.baseChange = baseAlgebra.inverse(Matrix.fromColumns(generators));
////		}
////
////	}
////
////	@Override
////	public LocalField<Base, Reduced> getBaseField() {
////		return base;
////	}
////
////	private Ext<Base> fromMultivariatePolynomialNoReduce(Polynomial<Base> e) {
////		return new Ext<>(base, unramified.degree(), totallyRamified.degree(), baseAlgebra, baseChange, e);
////	}
////
////	public Ext<Base> fromMultivariatePolynomial(Polynomial<Base> e) {
////		Monomial lc = e.leadingMonomial();
////		if (lc.exponents()[0] < unramified.degree() && lc.exponents()[1] < totallyRamified.degree()) {
////			return fromMultivariatePolynomialNoReduce(e);
////		}
////		if (lc.exponents()[0] < 2 * unramified.degree() && lc.exponents()[1] < 2 * totallyRamified.degree()) {
////			Polynomial<Base> reduced = ring.zero();
////			for (Monomial m : e.monomials()) {
////				reduced = ring.add(reduced, ring.multiply(e.coefficient(m), powers.get(m)));
////			}
////			return fromMultivariatePolynomialNoReduce(reduced);
////		}
////		return fromMultivariatePolynomialNoReduce(ideal.residue(e));
////	}
////
////	public PolynomialRing<Base> getPolynomialRing() {
////		return ring;
////	}
////
////	public Ext<Base> zero() {
////		return zero;
////	}
////
////	public Ext<Base> one() {
////		return one;
////	}
////
////	@Override
////	public Ext<Base> uniformizer() {
////		return uniformizer;
////	}
////
////	@Override
////	public Ext<Base> alpha() {
////		return add(reductionGenerator, uniformizer);
////	}
////
////	public Ext<Base> add(Ext<Base> t1, Ext<Base> t2) {
////		return fromMultivariatePolynomialNoReduce(ring.add(t1.e, t2.e));
////	}
////
////	public Ext<Base> negative(Ext<Base> t) {
////		return fromMultivariatePolynomialNoReduce(ring.negative(t.e));
////	}
////
////	public Ext<Base> multiply(Ext<Base> t1, Ext<Base> t2) {
////		return fromMultivariatePolynomial(ring.multiply(t1.e, t2.e));
////	}
////
////	public Ext<Base> inverse(Ext<Base> t) {
////		UnivariatePolynomialRing<Base> u = base.getUnivariatePolynomialRing();
////		UnivariatePolynomialRing<Polynomial<Base>> r = u.getUnivariatePolynomialRing();
////		UnivariatePolynomial<Polynomial<Base>> asUnivariate = ring.asUnivariatePolynomial(t.e, 2);
////		ExtendedResultantResult<Polynomial<Base>> ramifiedResult = r.extendedResultant(asUnivariate,
////				ring.asUnivariatePolynomial(totallyRamifiedMinimalPolynomial, 2));
////		if (ramifiedResult.getResultant().equals(u.zero())) {
////			throw new ArithmeticException("Division by zero!");
////		}
////		ExtendedResultantResult<Base> unramifiedResult = u.extendedResultant(ramifiedResult.getResultant(),
////				unramified.minimalPolynomial());
////		return fromMultivariatePolynomial(
////				ring.multiply(ring.getEmbedding(base.inverse(unramifiedResult.getResultant())),
////						ring.fromUnivariatePolynomial(ramifiedResult.getCoeff1(), 2),
////						ring.getEmbedding(unramifiedResult.getCoeff1(), new int[] { 0 })));
////	}
////
////	@Override
////	public FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> getEmbeddedExtension(
////			UnivariatePolynomial<Ext<Base>> minimalPolynomial) {
////		return LocalFieldExtension.createLocalFieldExtension(minimalPolynomial, this);
////	}
////	
////	@Override
////	public ReductionExtensionType reduction() {
////		return unramified.reduction();
////	}
////
//
//	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> createLocalFieldExtension(
//			UnivariatePolynomial<Base> minimalPolynomial, LocalField<Base, Reduced> base,
//			ReductionExtensionType reductionAsTrivialExtension) {
//		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> trivialExtension = new LocalFieldExtension<>(
//				base, reductionAsTrivialExtension);
//		return createLocalFieldExtension(trivialExtension.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial,
//				trivialExtension.getEmbeddingMap()), trivialExtension);
//	}
//
//	static <Base extends Element<Base>, Reduced extends Element<Reduced>, ReducedExt extends ExtensionFieldElement<Reduced, ReducedExt>, ReductionExtensionType extends FieldExtension<Reduced, ReducedExt, ReductionExtensionType>> LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> createLocalFieldExtension(
//			UnivariatePolynomial<Ext<Base>> minimalPolynomial,
//			LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> base) {
//		if (minimalPolynomial.degree() == 1) {
//			return base;
//		}
//		if (!base.isIrreducible(minimalPolynomial)) {
//			throw new ArithmeticException("Polynomial not irreducible!");
//		}
//		LocalRing<Ext<Base>, ReducedExt> baseRing = base.ringOfIntegers();
//		UnivariatePolynomialRing<ReducedExt> reducedPolynomials = base.reduction().getUnivariatePolynomialRing();
//		minimalPolynomial = baseRing.integralMinimalPolynomial(minimalPolynomial);
//		UnivariatePolynomial<ReducedExt> reduced = reducedPolynomials.toUnivariate(
//				reducedPolynomials.squareFreeFactorization(baseRing.reduceUnivariatePolynomial(minimalPolynomial))
//						.keySet().iterator().next());
//		if (reduced.degree() == minimalPolynomial.degree()) {
//			FieldEmbedding<Reduced, ReducedExt, ReductionExtensionType> reducedExtension = base.reduction()
//					.getEmbeddedExtension(reduced);
//			reduced = reducedExtension.minimalPolynomial();
//			return new LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>(
//					base.getBaseField().ringOfIntegers().liftUnivariatePolynomial(reduced), base.getBaseField(),
//					base.reduction);
//		}
//		if (baseRing.isEisenstein(minimalPolynomial)) {
//			return new LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>(minimalPolynomial, base);
//		}
//		LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> unramifiedExtension = base;
//		LocalRing<Ext<Base>, ReducedExt> unramifiedExtensionIntegers = base.ringOfIntegers();
//		UnivariatePolynomialRing<Ext<Base>> unramifiedPolynomials = base.getUnivariatePolynomialRing();
//		UnivariatePolynomial<Ext<Base>> minimalPolynomialOverUnramifiedExtension = minimalPolynomial;
//		while (true) {
//			if (minimalPolynomialOverUnramifiedExtension.degree() == 1) {
//				return unramifiedExtension;
//			} else if (reduced.degree() != 1) {
//				unramifiedExtension = unramifiedExtension
//						.getExtension(unramifiedExtensionIntegers.liftUnivariatePolynomial(reduced)).extension();
//				unramifiedExtensionIntegers = unramifiedExtension.ringOfIntegers();
//				unramifiedPolynomials = unramifiedExtension.getUnivariatePolynomialRing();
//				minimalPolynomialOverUnramifiedExtension = unramifiedPolynomials
//						.toUnivariate(unramifiedExtension
//								.factorization(unramifiedPolynomials.getEmbedding(minimalPolynomial,
//										unramifiedExtension.getEmbeddingMap()))
//								.getFactors().keySet().iterator().next());
//				reducedPolynomials = unramifiedExtension.reduction().getUnivariatePolynomialRing();
//				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
//						.squareFreeFactorization(unramifiedExtensionIntegers
//								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
//						.keySet().iterator().next());
//			} else if (!reduced.univariateCoefficient(0).equals(unramifiedExtension.reduction().zero())) {
//				Ext<Base> constant = unramifiedExtension.lift(reduced.univariateCoefficient(0));
//				Polynomial<Ext<Base>> transform = unramifiedPolynomials
//						.getPolynomial(unramifiedExtension.negative(constant), unramifiedExtension.one());
//				minimalPolynomialOverUnramifiedExtension = unramifiedPolynomials
//						.substitute(minimalPolynomialOverUnramifiedExtension, Collections.singletonList(transform));
//				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
//						.squareFreeFactorization(unramifiedExtensionIntegers
//								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
//						.keySet().iterator().next());
//			} else if (unramifiedExtension.valuation(minimalPolynomialOverUnramifiedExtension.univariateCoefficient(0))
//					.value() == 1) {
//				FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> totallyRamifiedExtension = unramifiedExtension
//						.getEmbeddedExtension(minimalPolynomialOverUnramifiedExtension);
//				Ext<Base> alpha = totallyRamifiedExtension.getField().add(
//						((LocalField<Ext<Base>, ReducedExt>) totallyRamifiedExtension.getField()).uniformizer(),
//						totallyRamifiedExtension.getEmbedding(unramifiedExtension.alpha()));
//				minimalPolynomial = totallyRamifiedExtension.minimalPolynomial(alpha);
//				return new LocalFieldExtension<>(minimalPolynomial, base);
//			} else {
//				minimalPolynomialOverUnramifiedExtension = unramifiedExtensionIntegers
//						.integralMinimalPolynomial(minimalPolynomialOverUnramifiedExtension);
//				reduced = reducedPolynomials.toUnivariate(reducedPolynomials
//						.squareFreeFactorization(unramifiedExtensionIntegers
//								.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
//						.keySet().iterator().next());
//				if (unramifiedExtension.valuation(minimalPolynomialOverUnramifiedExtension.univariateCoefficient(0))
//						.value() > 1) {
//					GenericExtensionField<Ext<Base>> extension = new GenericExtensionField<>(
//							minimalPolynomialOverUnramifiedExtension, unramifiedExtension);
//					int ramificationIndex = minimalPolynomialOverUnramifiedExtension.degree() / reduced.degree();
//					GenericExtensionFieldElement<Ext<Base>> alpha = extension.alpha();
//					int alphaValuation = unramifiedExtension.valuation(extension.norm(alpha)).value()
//							/ reduced.degree();
//					while (true) {
//						if (BigInteger.valueOf(alphaValuation).gcd(BigInteger.valueOf(ramificationIndex))
//								.intValueExact() == 1) {
//							int power = BigInteger.valueOf(alphaValuation)
//									.modInverse(BigInteger.valueOf(ramificationIndex)).intValueExact();
//							int ramificationPower = (power * alphaValuation - 1) / ramificationIndex;
//							minimalPolynomialOverUnramifiedExtension = extension
//									.minimalPolynomial(extension.divide(extension.power(alpha, power),
//											extension.power(extension.getEmbedding(unramifiedExtension.uniformizer()),
//													ramificationPower)));// ,
//							reduced = reducedPolynomials.toUnivariate(reducedPolynomials
//									.squareFreeFactorization(unramifiedExtensionIntegers
//											.reduceUnivariatePolynomial(minimalPolynomialOverUnramifiedExtension))
//									.keySet().iterator().next());
//							break;
//						} else {
//							int gcd = BigInteger.valueOf(alphaValuation).gcd(BigInteger.valueOf(ramificationIndex))
//									.intValueExact();
//							alpha = extension.add(extension.power(alpha, ramificationIndex / gcd),
//									extension.getEmbedding(unramifiedExtension.power(unramifiedExtension.uniformizer(),
//											alphaValuation / gcd)));
//							alphaValuation = unramifiedExtension.valuation(extension.norm(alpha)).value()
//									/ reduced.degree();
//						}
//					}
//				}
//			}
//		}
//	}
//
//	private LocalFieldExtension(UnivariatePolynomial<Base> minimalPolynomial, LocalField<Base, Reduced> base,
//			ReductionExtensionType reductionAsTrivialExtension) {
//		super(minimalPolynomial, base);
//		// this.base = base;
//		this.prime = base;
//		this.highAccuracyPrime = prime.withAccuracy(prime.getAccuracy() * degree());
//		this.highAccuracyPolynomials = highAccuracyPrime.getUnivariatePolynomialRing();
//		this.highAccuracyAlgebra = new FreeModule<>(highAccuracyPrime, degree()).matrixAlgebra();
//		this.conjugates = new TreeMap<>();
//		this.extensions = new TreeMap<>();
//		this.localRing = new LocalRingImplementation<>(this,
//				base.ringOfIntegers().toString() + "[X]/(" + minimalPolynomial + ")");
//		// this.trivial = false;
//		LocalRing<Base, Reduced> baseRing = base.ringOfIntegers();
//		UnivariatePolynomialRing<Reduced> reducedPolynomials = base.reduction().getUnivariatePolynomialRing();
//		UnivariatePolynomial<Reduced> reduced = reducedPolynomials.toUnivariate(
//				reducedPolynomials.squareFreeFactorization(baseRing.reduceUnivariatePolynomial(minimalPolynomial))
//						.keySet().iterator().next());
//		this.reduction = reductionAsTrivialExtension.getExtension(reductionAsTrivialExtension
//				.getUnivariatePolynomialRing().getEmbedding(reduced, reductionAsTrivialExtension.getEmbeddingMap()))
//				.extension();
//		this.residueDegree = reduction().degree();
//		this.ramificationIndex = degree() / residueDegree;
//		UnivariatePolynomial<Ext<Base>> residuePolynomial = getUnivariatePolynomialRing()
//				.getEmbedding(baseRing.liftUnivariatePolynomial(reduced), getEmbeddingMap());
//		this.uniformizer = getUnivariatePolynomialRing().evaluate(residuePolynomial, alpha());
//		if (this.uniformizer.equals(zero())) {
//			this.uniformizer = getEmbedding(base.uniformizer());
//		}
//		UnivariatePolynomial<Ext<Base>> lifted = getUnivariatePolynomialRing().getEmbedding(
//				prime.ringOfIntegers().liftUnivariatePolynomial(reduction.minimalPolynomial()), getEmbeddingMap());
//		this.liftOfReductionGenerator = localRing.henselLift(lifted, reduction.alpha());
//		List<Vector<Base>> reductionGeneratorBase = new ArrayList<>();
//		reductionGeneratorUniformizerBase = new ArrayList<>();
//		for (int i = 0; i < ramificationIndex; i++) {
//			for (int j = 0; j < residueDegree; j++) {
//				Ext<Base> generator = multiply(power(liftOfReductionGenerator, j), power(uniformizer, i));
//				reductionGeneratorUniformizerBase.add(generator);
//				reductionGeneratorBase.add(asVector(generator));
//			}
//		}
//		this.gammaToReductionUniformizer = matrixAlgebra().inverse(Matrix.fromColumns(reductionGeneratorBase));
//		if (prime.valuation(norm(uniformizer)).value() != residueDegree) {
//			throw new ArithmeticException("preparation wrong (uniformizer is not valuation 1 over prime)");
//		}
//	}
//
//	@Override
//	public Exactness exactness() {
//		return prime.exactness();
//	}
//
//	@Override
//	public FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> getEmbeddedExtension(
//			UnivariatePolynomial<Ext<Base>> minimalPolynomial) {
//
//		FieldEmbedding<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> extension;
//		if (extensions.containsKey(minimalPolynomial)) {
//			extension = extensions.get(minimalPolynomial);
//		} else {
//			extension = createLocalFieldExtension(minimalPolynomial, this);
//			extensions.put(minimalPolynomial, extension);
//		}
//		return extension;
//	}
//
//	public LocalField<Base, Reduced> getBaseField() {
//		return prime;
//	}
//
//	public int ramificationIndex() {
//		return ramificationIndex;
//	}
//
//	public int residueDegree() {
//		return residueDegree;
//	}
//
//	@Override
//	public Ext<Base> negative(Ext<Base> t, int accuracy) {
//		if (accuracy < this.getAccuracy()) {
//			return round(negative(t), accuracy);
//		} else if (accuracy == this.getAccuracy()) {
//			return negative(t);
//		}
//		return withAccuracy(accuracy).negative(t);
//	}
//
//	@Override
//	public Ext<Base> inverse(Ext<Base> t) {
//		UnivariatePolynomialRing.ExtendedResultantResult<Base> eres = highAccuracyPolynomials
//				.extendedResultant(minimalPolynomial(), t.e);
//		if (eres.getResultant().equals(prime.zero())) {
//			throw new ArithmeticException("Not invertible!");
//		}
//		return fromPolynomial(prime.ringOfIntegers().roundUnivariatePolynomial(
//				highAccuracyPolynomials.divideScalar(eres.getCoeff2(), eres.getGcd().leadingCoefficient()),
//				prime.getAccuracy()));
//	}
//
//	@Override
//	public Ext<Base> inverse(Ext<Base> t, int accuracy) {
//		if (accuracy < this.getAccuracy()) {
//			return round(inverse(t), accuracy);
//		} else if (accuracy == this.getAccuracy()) {
//			return inverse(t);
//		}
//		return withAccuracy(accuracy).inverse(t);
//	}
//
//	@Override
//	public Ext<Base> divide(Ext<Base> dividend, Ext<Base> divisor) {
//		if (dividend.equals(zero())) {
//			return zero();
//		}
//		if (divisor.equals(zero())) {
//			throw new ArithmeticException("Division by zero!");
//		}
//		return multiply(dividend, inverse(divisor,
//				getAccuracy() + Math.max(0, -valuation(dividend).value() + valuation(divisor).value())));
//	}
//
//	/*
//	 * public Ext<Base, Reduced, ReducedExt> fromRational(Fraction t) { return
//	 * getPrimeEmbedding(prime.fromRational(t)); }
//	 */
//	@Override
//	public Ext<Base> round(Ext<Base> t, int accuracy) {
//		Vector<Base> asVector = matrixAlgebra().multiply(gammaToReductionUniformizer, asVector(t));
//		Ext<Base> rounded = zero();
//		for (int i = 0; i < ramificationIndex; i++) {
//			for (int j = 0; j < residueDegree; j++) {
//				Ext<Base> generator = reductionGeneratorUniformizerBase.get(i * residueDegree + j);
//				rounded = add(rounded, multiply(getEmbedding(prime.round(asVector.get(residueDegree * i + j + 1),
//						(accuracy - i + ramificationIndex - 1) / ramificationIndex)), generator));
//			}
//		}
//		return rounded;
//	}
//
//	public ReductionExtensionType reduction() {
//		return reduction;
//	}
//
//	public ReducedExt reduce(Ext<Base> t) {
//		if (liftOfReductionGenerator != null) {
//			Vector<Base> asVector = asVector(t);
//			Vector<Base> asVectorReductionUniformizer = matrixAlgebra().multiply(gammaToReductionUniformizer, asVector);
//			Polynomial<Base> asPolynomial = prime.getUnivariatePolynomialRing()
//					.getPolynomial(asVectorReductionUniformizer.asList().subList(0, residueDegree));
//			UnivariatePolynomial<Reduced> asReducedPolynomial = prime.reduction().getUnivariatePolynomialRing()
//					.getEmbedding(asPolynomial, new MathMap<>() {
//						@Override
//						public Reduced evaluate(Base t) {
//							return prime.reduce(t);
//						}
//					});
//			UnivariatePolynomialRing<ReducedExt> reducedRing = reduction().getUnivariatePolynomialRing();
//			return reducedRing.evaluate(reducedRing.getEmbedding(asReducedPolynomial, reduction().getEmbeddingMap()),
//					reduction().alpha());
//		}
//		UnivariatePolynomial<Base> asPolynomial = t.asPolynomial();
//		UnivariatePolynomialRing<ReducedExt> reducedRing = reduction().getUnivariatePolynomialRing();
//		return reducedRing
//				.evaluate(reducedRing.getEmbedding(prime.ringOfIntegers().reduceUnivariatePolynomial(asPolynomial),
//						reduction().getEmbeddingMap()), reduction().alpha());
//	}
//
//	public Ext<Base> lift(ReducedExt t) {
//		UnivariatePolynomialRing<Ext<Base>> ring = getUnivariatePolynomialRing();
//		UnivariatePolynomial<Reduced> asReducedPolynomial = t.asPolynomial();
//		UnivariatePolynomial<Base> asPolynomial = prime.getUnivariatePolynomialRing().getEmbedding(asReducedPolynomial,
//				new MathMap<>() {
//					@Override
//					public Base evaluate(Reduced t) {
//						return prime.lift(t);
//					}
//				});
//		if (liftOfReductionGenerator != null) {
//			return ring.evaluate(ring.getEmbedding(asPolynomial, getEmbeddingMap()), liftOfReductionGenerator);
//		}
//		return ring.evaluate(ring.getEmbedding(asPolynomial, getEmbeddingMap()), alpha());
//
//	}
//
//	@Override
//	public boolean isIrreducible(UnivariatePolynomial<Ext<Base>> t) {
//		return ringOfIntegers().isIrreducible(t);
//	}
//
//	@Override
//	public FactorizationResult<Polynomial<Ext<Base>>> factorization(UnivariatePolynomial<Ext<Base>> t) {
//		return ringOfIntegers().factorization(t, true);
//	}
//
//	@Override
//	public Ext<Base> getRandomInteger() {
//		List<Base> coefficients = new ArrayList<>();
//		for (int i = 0; i < degree(); i++) {
//			coefficients.add(prime.getRandomInteger());
//		}
//		return fromPolynomial(prime.getUnivariatePolynomialRing().getPolynomial(coefficients));
//	}
//
//	@Override
//	public int getAccuracy() {
//		return prime.getAccuracy() * ramificationIndex;
//	}
//
//	@Override
//	public LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> withAccuracy(int accuracy) {
//		if (this.degree() == 1) {
//			return new LocalFieldExtension<>(prime.withAccuracy(accuracy), reduction);
//		}
//		LocalField<Base, Reduced> base = this.prime.withAccuracy(accuracy);
//		UnivariatePolynomial<Base> minimalPolynomial = base.ringOfIntegers()
//				.roundUnivariatePolynomial(minimalPolynomial(), accuracy);
//		return new LocalFieldExtension<>(minimalPolynomial, base, reduction);
//	}
//
//	@Override
//	public GaloisGroup<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>> galoisGroup() {
//		if (ramificationIndex == 1) {
//			List<FieldAutomorphism<Base, Ext<Base>, LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType>>> elements = new ArrayList<>();
//			for (FieldAutomorphism<Reduced, ReducedExt, ReductionExtensionType> reduced : reduction().galoisGroup()) {
//				elements.add(new FieldAutomorphism<>(this, reduced.asArray()));
//			}
//			return new GaloisGroup<>(this, elements);
//		}
//		throw new UnsupportedOperationException("Not implemented");
//	}
//
//	@Override
//	public List<Ext<Base>> conjugates(Ext<Base> s) {
//		if (conjugates.containsKey(s)) {
//			return conjugates.get(s);
//		}
//		UnivariatePolynomial<Ext<Base>> mipo = getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(s),
//				getEmbeddingMap());
//		if (ramificationIndex == 1) {
//			List<Ext<Base>> conjugates = new ArrayList<>();
//			for (ReducedExt c : reduction().conjugates(reduce(s))) {
//				conjugates.add(localRing.henselLift(mipo, c));
//			}
//			this.conjugates.put(s, conjugates);
//			return conjugates;
//		}
//		List<Ext<Base>> conjugates = new ArrayList<>();
//		conjugates.addAll(roots(mipo).keySet());
//		this.conjugates.put(s, conjugates);
//		return conjugates(s);
//	}
//
//	@Override
//	public boolean isInteger(Ext<Base> t) {
//		return valuation(t).compareTo(new Value(0)) >= 0;
//	}
//
//	@Override
//	public double value(Ext<Base> t) {
//		return Math.pow(prime.value(norm(t)), 1.0 / degree());
//	}
//
//	@Override
//	public Value valuation(Ext<Base> t) {
//		if (t.valuation != null) {
//			return t.valuation;
//		}
//		Base norm = highAccuracyAlgebra.determinant(asMatrix(t));
////		if (!t.equals(zero()) && normOverPrime.equals(prime.zero())) {
////			int minValuation = -1;
////			for (int i = 0; i < totalDegree(); i++) {
////				Base c = t.e.asPolynomial().univariateCoefficient(i);
////				if (c.equals(prime.zero())) {
////					continue;
////				}
////				int valuation = prime.valuation(c);
////				if (minValuation == -1 || valuation < minValuation) {
////					minValuation = valuation;
////				}
////			}
////			return valuation(multiply(t, power(inverse(uniformizer()), minValuation))) + minValuation;
////		}
//		Value normValue = highAccuracyPrime.valuation(norm);
//		t.valuation = normValue.isInfinite() ? normValue : new Value(normValue.value() / residueDegree);
//		return t.valuation;
//	}
//
//	@Override
//	public Ext<Base> uniformizer() {
//		return uniformizer;
//	}
//
//	@Override
//	public LocalRing<Ext<Base>, ReducedExt> ringOfIntegers() {
//		return localRing;
//	}
//
//	@Override
//	protected Ext<Base> fromSmallDegreePolynomial(UnivariatePolynomial<Base> polynomial) {
//		return new Ext<Base>(polynomial);
//	}
//
//	@Override
//	public LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> makeExtension(
//			UnivariatePolynomial<Base> minimalPolynomial) {
//		return createLocalFieldExtension(minimalPolynomial, prime, reduction);
//	}
//
//	@Override
//	protected LocalFieldExtension<Base, Reduced, ReducedExt, ReductionExtensionType> asExtensionType() {
//		return this;
//	}
//}
