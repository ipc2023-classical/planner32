(define (domain sat3-dom-65-15-1390490377005624047)
  (:requirements :typing :strips)
  (:types var clause)

  (:constants
    v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 - var
    c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33 c34 c35 c36 c37 c38 c39 c40 c41 c42 c43 c44 c45 c46 c47 c48 c49 c50 c51 c52 c53 c54 c55 c56 c57 c58 c59 c60 c61 c62 c63 c64 - clause
  )

  (:predicates
    (solved ?c - clause)
    (unassigned ?v - var)
  )

  (:action assign-v1-true
    :parameters ()
    :precondition (unassigned v1)
    :effect
    (and
      (not (unassigned v1))
      (solved c23)
      (solved c35)
      (solved c63)
      (solved c64)
    )
  )

  (:action assign-v2-true
    :parameters ()
    :precondition (unassigned v2)
    :effect
    (and
      (not (unassigned v2))
      (solved c2)
      (solved c9)
      (solved c12)
      (solved c15)
      (solved c20)
      (solved c21)
      (solved c28)
      (solved c43)
      (solved c45)
      (solved c52)
    )
  )

  (:action assign-v3-true
    :parameters ()
    :precondition (unassigned v3)
    :effect
    (and
      (not (unassigned v3))
      (solved c2)
      (solved c18)
      (solved c21)
      (solved c36)
      (solved c40)
      (solved c42)
      (solved c46)
      (solved c62)
    )
  )

  (:action assign-v4-true
    :parameters ()
    :precondition (unassigned v4)
    :effect
    (and
      (not (unassigned v4))
      (solved c13)
      (solved c14)
      (solved c35)
      (solved c39)
      (solved c44)
      (solved c50)
      (solved c53)
    )
  )

  (:action assign-v5-true
    :parameters ()
    :precondition (unassigned v5)
    :effect
    (and
      (not (unassigned v5))
      (solved c25)
      (solved c51)
      (solved c61)
      (solved c63)
    )
  )

  (:action assign-v6-true
    :parameters ()
    :precondition (unassigned v6)
    :effect
    (and
      (not (unassigned v6))
      (solved c8)
      (solved c20)
      (solved c24)
      (solved c34)
      (solved c49)
      (solved c59)
    )
  )

  (:action assign-v7-true
    :parameters ()
    :precondition (unassigned v7)
    :effect
    (and
      (not (unassigned v7))
      (solved c2)
      (solved c10)
      (solved c16)
      (solved c26)
      (solved c64)
    )
  )

  (:action assign-v8-true
    :parameters ()
    :precondition (unassigned v8)
    :effect
    (and
      (not (unassigned v8))
      (solved c16)
      (solved c18)
      (solved c19)
      (solved c30)
      (solved c33)
      (solved c37)
      (solved c43)
      (solved c52)
      (solved c59)
      (solved c61)
    )
  )

  (:action assign-v9-true
    :parameters ()
    :precondition (unassigned v9)
    :effect
    (and
      (not (unassigned v9))
      (solved c10)
      (solved c25)
      (solved c34)
      (solved c35)
      (solved c37)
      (solved c39)
      (solved c46)
      (solved c47)
      (solved c51)
      (solved c53)
      (solved c57)
    )
  )

  (:action assign-v10-true
    :parameters ()
    :precondition (unassigned v10)
    :effect
    (and
      (not (unassigned v10))
      (solved c11)
      (solved c21)
      (solved c27)
      (solved c31)
      (solved c33)
      (solved c40)
    )
  )

  (:action assign-v11-true
    :parameters ()
    :precondition (unassigned v11)
    :effect
    (and
      (not (unassigned v11))
      (solved c1)
      (solved c11)
      (solved c56)
      (solved c60)
    )
  )

  (:action assign-v12-true
    :parameters ()
    :precondition (unassigned v12)
    :effect
    (and
      (not (unassigned v12))
      (solved c4)
      (solved c6)
      (solved c19)
      (solved c41)
      (solved c46)
      (solved c53)
    )
  )

  (:action assign-v13-true
    :parameters ()
    :precondition (unassigned v13)
    :effect
    (and
      (not (unassigned v13))
      (solved c22)
      (solved c56)
      (solved c57)
      (solved c59)
      (solved c63)
    )
  )

  (:action assign-v14-true
    :parameters ()
    :precondition (unassigned v14)
    :effect
    (and
      (not (unassigned v14))
      (solved c5)
      (solved c6)
      (solved c8)
      (solved c20)
      (solved c32)
      (solved c43)
      (solved c55)
    )
  )

  (:action assign-v15-true
    :parameters ()
    :precondition (unassigned v15)
    :effect
    (and
      (not (unassigned v15))
      (solved c0)
      (solved c1)
      (solved c22)
      (solved c24)
      (solved c26)
      (solved c31)
      (solved c32)
    )
  )

  (:action assign-v1-false
    :parameters ()
    :precondition (unassigned v1)
    :effect
    (and
      (not (unassigned v1))
      (solved c1)
      (solved c4)
      (solved c8)
      (solved c47)
      (solved c55)
    )
  )

  (:action assign-v2-false
    :parameters ()
    :precondition (unassigned v2)
    :effect
    (and
      (not (unassigned v2))
      (solved c0)
      (solved c23)
      (solved c24)
      (solved c27)
      (solved c40)
      (solved c41)
      (solved c55)
      (solved c58)
      (solved c60)
    )
  )

  (:action assign-v3-false
    :parameters ()
    :precondition (unassigned v3)
    :effect
    (and
      (not (unassigned v3))
      (solved c5)
      (solved c9)
      (solved c57)
      (solved c58)
      (solved c61)
    )
  )

  (:action assign-v4-false
    :parameters ()
    :precondition (unassigned v4)
    :effect
    (and
      (not (unassigned v4))
      (solved c4)
      (solved c5)
      (solved c7)
      (solved c12)
      (solved c29)
      (solved c49)
      (solved c52)
    )
  )

  (:action assign-v5-false
    :parameters ()
    :precondition (unassigned v5)
    :effect
    (and
      (not (unassigned v5))
      (solved c7)
      (solved c9)
      (solved c15)
      (solved c32)
      (solved c34)
      (solved c36)
      (solved c58)
      (solved c62)
    )
  )

  (:action assign-v6-false
    :parameters ()
    :precondition (unassigned v6)
    :effect
    (and
      (not (unassigned v6))
      (solved c6)
      (solved c10)
      (solved c44)
      (solved c48)
      (solved c64)
    )
  )

  (:action assign-v7-false
    :parameters ()
    :precondition (unassigned v7)
    :effect
    (and
      (not (unassigned v7))
      (solved c13)
      (solved c14)
      (solved c15)
      (solved c19)
      (solved c42)
      (solved c47)
      (solved c48)
      (solved c49)
      (solved c54)
      (solved c60)
    )
  )

  (:action assign-v8-false
    :parameters ()
    :precondition (unassigned v8)
    :effect
    (and
      (not (unassigned v8))
      (solved c14)
      (solved c28)
      (solved c29)
      (solved c31)
      (solved c50)
    )
  )

  (:action assign-v9-false
    :parameters ()
    :precondition (unassigned v9)
    :effect
    (and
      (not (unassigned v9))
      (solved c11)
      (solved c17)
      (solved c26)
      (solved c27)
      (solved c28)
      (solved c29)
      (solved c45)
      (solved c48)
      (solved c50)
      (solved c62)
    )
  )

  (:action assign-v10-false
    :parameters ()
    :precondition (unassigned v10)
    :effect
    (and
      (not (unassigned v10))
      (solved c3)
      (solved c13)
      (solved c23)
      (solved c30)
      (solved c38)
    )
  )

  (:action assign-v11-false
    :parameters ()
    :precondition (unassigned v11)
    :effect
    (and
      (not (unassigned v11))
      (solved c36)
      (solved c38)
      (solved c42)
      (solved c45)
    )
  )

  (:action assign-v12-false
    :parameters ()
    :precondition (unassigned v12)
    :effect
    (and
      (not (unassigned v12))
      (solved c0)
      (solved c17)
      (solved c37)
      (solved c39)
      (solved c54)
    )
  )

  (:action assign-v13-false
    :parameters ()
    :precondition (unassigned v13)
    :effect
    (and
      (not (unassigned v13))
      (solved c3)
      (solved c16)
      (solved c18)
      (solved c44)
      (solved c54)
    )
  )

  (:action assign-v14-false
    :parameters ()
    :precondition (unassigned v14)
    :effect
    (and
      (not (unassigned v14))
      (solved c3)
      (solved c7)
      (solved c12)
      (solved c17)
      (solved c22)
      (solved c30)
      (solved c33)
      (solved c51)
    )
  )

  (:action assign-v15-false
    :parameters ()
    :precondition (unassigned v15)
    :effect
    (and
      (not (unassigned v15))
      (solved c25)
      (solved c38)
      (solved c41)
      (solved c56)
    )
  )

)
