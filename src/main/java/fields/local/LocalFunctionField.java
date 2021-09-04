//package fields.local;
//
//import java.math.BigInteger;
//import java.util.Iterator;
//
//import fields.exceptions.InfinityException;
//import fields.floatingpoint.Reals.Real;
//import fields.helper.AbstractField;
//import fields.helper.CoordinateRing.CoordinateRingElement;
//import fields.interfaces.Element;
//import fields.interfaces.Field;
//import fields.interfaces.Ideal;
//import fields.interfaces.LocalField;
//import fields.interfaces.LocalRing;
//import fields.local.FormalPowerSeries.PowerSeries;
//import fields.vectors.pivot.PivotStrategy;
//import fields.vectors.pivot.ValuationPivotStrategy;
//import util.Identity;
//import varieties.FunctionField;
//import varieties.RationalFunctionTmp;
//
//public class LocalFunctionField<T extends Element<T>> extends AbstractField<RationalFunctionTmp<T>>
//		implements LocalField<RationalFunctionTmp<T>, T> {
//	private FunctionField<T> field;
//	private Ideal<CoordinateRingElement<T>> maximalIdeal;
//	private LocalRing<RationalFunctionTmp<T>, T> localRing;
//
//	public LocalFunctionField(FunctionField<T> field, Ideal<CoordinateRingElement<T>> maximalIdeal) {
//		this.field = field;
//		this.maximalIdeal = maximalIdeal;
//		this.localRing = new LocalRingImplementation<>(this, field.getCoordinateRing().toString());
//	}
//
//	@Override
//	public Real value(RationalFunctionTmp<T> t) {
//		return null;
//	}
//
//	@Override
//	public RationalFunctionTmp<T> zero() {
//		return field.zero();
//	}
//
//	@Override
//	public RationalFunctionTmp<T> one() {
//		return field.one();
//	}
//
//	@Override
//	public BigInteger characteristic() {
//		return field.characteristic();
//	}
//
//	@Override
//	public RationalFunctionTmp<T> add(RationalFunctionTmp<T> t1, RationalFunctionTmp<T> t2) {
//		return field.add(t1, t2);
//	}
//
//	@Override
//	public RationalFunctionTmp<T> negative(RationalFunctionTmp<T> t) {
//		return field.negative(t);
//	}
//
//	@Override
//	public RationalFunctionTmp<T> multiply(RationalFunctionTmp<T> t1, RationalFunctionTmp<T> t2) {
//		return field.multiply(t1, t2);
//	}
//
//	@Override
//	public RationalFunctionTmp<T> inverse(RationalFunctionTmp<T> t) {
//		return field.inverse(t);
//	}
//
//	@Override
//	public Exactness exactness() {
//		return field.exactness();
//	}
//
//	@Override
//	public RationalFunctionTmp<T> getRandomElement() {
//		return field.getRandomElement();
//	}
//
//	@Override
//	public boolean isFinite() {
//		return false;
//	}
//
//	@Override
//	public BigInteger getNumberOfElements() throws InfinityException {
//		throw new InfinityException();
//	}
//
//	@Override
//	public Iterator<RationalFunctionTmp<T>> iterator() {
//		return field.iterator();
//	}
//
//	@Override
//	public boolean isComplete() {
//		return false;
//	}
//
//	@Override
//	public RationalFunctionTmp<T> inverse(RationalFunctionTmp<T> t, int accuracy) {
//		return inverse(t);
//	}
//	
//	@Override
//	public PivotStrategy<RationalFunctionTmp<T>> preferredPivotStrategy() {
//		return new ValuationPivotStrategy<>(this.valuation());
//	}
//
//	@Override
//	public RationalFunctionTmp<T> negative(RationalFunctionTmp<T> t, int accuracy) {
//		return negative(t);
//	}
//
//	@Override
//	public boolean isInteger(RationalFunctionTmp<T> t) {
//		// TODO Auto-generated method stub
//		return false;
//	}
//
//	@Override
//	public Value valuation(RationalFunctionTmp<T> t) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public RationalFunctionTmp<T> uniformizer() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public Field<T> reduction() {
//		return field.getField();
//	}
//
//	@Override
//	public T reduce(RationalFunctionTmp<T> t) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//	
//	@Override
//	public RationalFunctionTmp<T> upToUniformizerPower(RationalFunctionTmp<T> t) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public RationalFunctionTmp<T> lift(T s) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public RationalFunctionTmp<T> round(RationalFunctionTmp<T> t, int accuracy) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public int getAccuracy() {
//	throw new InfinityException();}
//
//	@Override
//	public LocalField<RationalFunctionTmp<T>, T> withAccuracy(int accuracy) {
//		return this;
//	}
//
//	@Override
//	public OtherVersion<RationalFunctionTmp<T>,RationalFunctionTmp<T>, T> exact() {
//		return new OtherVersion<>(this, new Identity<>(), new Identity<>());
//	}
//
//	@Override
//	public OtherVersion<RationalFunctionTmp<T>,  PowerSeries<T>, T> complete(int accuracy) {
//		FormalPowerSeries<T> formalPowerSeries = new FormalPowerSeries<>(reduction(), accuracy);
//		return new OtherVersion<RationalFunctionTmp<T>,  PowerSeries<T>, T>(formalPowerSeries, formalPowerSeries.fromRationalFunctionMap(), formalPowerSeries.toRationalFunctionMap());
//	}
//
//	@Override
//	public RationalFunctionTmp<T> getRandomInteger() {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public LocalRing<RationalFunctionTmp<T>, T> ringOfIntegers() {
//		return localRing;
//	}
//
//}
