package fields.quaternions;

import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.quaternions.AbstractQuaternions.Quaternion;

public class RationalQuaternionOrder extends AbstractQuaternionOrder<Fraction, IntE, PFE> {

	public RationalQuaternionOrder(Quaternions<Fraction> quaternions, Quaternion<Fraction> t1, Quaternion<Fraction> t2,
			Quaternion<Fraction> t3) {
		super(Integers.z(), quaternions, t1, t2, t3);
	}
	
	
}
