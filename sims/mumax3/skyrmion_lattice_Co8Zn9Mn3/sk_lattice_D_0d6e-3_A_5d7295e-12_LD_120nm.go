//+build ignore

package main

import (
	. "github.com/mumax/3/engine"
    "math"
    "fmt"
)

func main() {

	defer InitAndClose()()

    cells := 5e-9
    Lx := 1500e-9
    Ly := 700e-9
    Lz := 200e-9

    lx := int(Lx / cells)
    ly := int(Ly / cells)
    lz := int(Lz / cells)
	SetGridSize(lx, ly, lz)
	SetCellSize(cells, cells, cells)

    Ac := 5.72957795e-12
    Dc := 0.6e-3
    Ms := 4.60545e5
    Msat.Set(Ms)
    Aex.Set(Ac)
    Dbulk.Set(Dc)
    // Alpha.Set(0.02)
    LD := 4 * math.Pi * (Ac) / (Dc)
    // BD := 0.5 * (Dc * Dc) / (Ac * Ms)

    // Conical phase angle params:
    // BDratio := 0.0 / BD

    // Magnetisation ---------------------------------------------------------
    
    L_sk := LD
    F := 2 * math.Pi / L_sk

    q0 := Vector(-1.0, 0.0, 0.0)
    q0 = q0.Mul(F)
    q1 := Vector(0.5, -math.Sqrt(3) * 0.5, 0.0)
    q1 = q1.Mul(F)
    q2 := Vector(0.5, math.Sqrt(3) * 0.5, 0.0)
    q2 = q2.Mul(F)

    e0 := Vector(0.0, 1.0, 0.0)
    e1 := Vector(-math.Sqrt(3) * 0.5, -0.5, 0.0)
    e2 := Vector(math.Sqrt(3) * 0.5, -0.5, 0.0)

    ez := Vector(0.0, 0.0, -1.0)

    host := M.Buffer().HostCopy()
    h := host.Vectors()
    n := M.Mesh().Size()

    // Set a function to update the magnetisation (using h and host array)
    // using the variables from the main function
    magFunction := func() {
        
        for iz := 0; iz < n[Z]; iz++ {
            for iy := 0; iy < n[Y]; iy++ {
                for ix := 0; ix < n[X]; ix++ {
                    pos := Index2Coord(ix, iy, iz)
                    x, y, z := pos[X], pos[Y], pos[Z]
                    r := Vector(x, y, z)

                    m := Vector(0.0, 0.0, 0.0)

                    for m_idx := 0; m_idx < 3; m_idx++ {
                        m[m_idx] = (ez[m_idx] * math.Cos(q0.Dot(r) + math.Pi) -
                                    e0[m_idx] * math.Sin(q0.Dot(r) + math.Pi) +
                                    ez[m_idx] * math.Cos(q1.Dot(r) + math.Pi) -
                                    e1[m_idx] * math.Sin(q1.Dot(r) + math.Pi) +
                                    ez[m_idx] * math.Cos(q2.Dot(r) + math.Pi) -
                                    e2[m_idx] * math.Sin(q2.Dot(r) + math.Pi))
                   }

                    h[X][iz][iy][ix] = float32(m[X])
                    h[Y][iz][iy][ix] = float32(m[Y])
                    h[Z][iz][iy][ix] = float32(m[Z])
                }
            }
        }
    }

    magFunction()
    M.SetArray(host)

    // -----------------------------------------------------------------------

    TableAdd(E_total)
	TableAdd(B_ext)
	TableAdd(E_exch)
	TableAdd(MaxTorque)

    // Relax with conjugate gradient:
    // MinimizerStop = 1e-6
    StopMaxDm = 1e-6
    // Relax with LLG:
    Alpha.Set(0.9)
    Precess = false
    RelaxTorqueThreshold = 1e-3
	// AutoSave(&M, 5e-10)
	TableAutoSave(5e-10)

    Bmin := -281e-3
    Bstep := 20e-3
    simname := fmt.Sprintf("")

    for B := -20e-3; B >= Bmin; B -= Bstep {
        
        // Reset magnetisation:
        magFunction()
        M.SetArray(host)
        // SaveAs(&M, "initial_state")

        // simname = fmt.Sprintf("m_skX_Bz_%06.0f_mT", B * 1000)
        // SaveAs(&M, simname)

        B_ext.Set(Vector(0, 0, B))

        // Relax with LLG and then minimize with SD
        Run(3e-9)
        TableSave()
        Minimize()   // small changes best minimized by minimize()
        TableSave()

        // Save the magnetisation:
        simname = fmt.Sprintf("m_skX_Bz_%06.0f_mT", B * 1000)
        SaveAs(&M, simname)
    }
}
