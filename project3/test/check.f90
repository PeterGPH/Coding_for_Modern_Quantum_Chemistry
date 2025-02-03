program test_mod_scf
    use mod_scf
    use iso_fortran_env, only: dp => real64
    implicit none

    type(scf) :: mol
    integer :: i, j, k, l, io_status
    character(len=100) :: filename1, filename2, filename3, filename4, filename5, filename6
    real(dp), parameter :: tolerance = 1.0e-4_dp
    ! Write test file
    filename1 = "geom.dat"
    filename2 = 'enuc.dat'
    filename3 = 's.dat'
    filename4 = 't.dat'
    filename5 = 'v.dat'
    filename6 = 'eri.dat'
    call write_test_file(filename1, filename2, filename3, filename4, filename5, filename6)

    ! Read geometry from file
    call mol%read_geometry(filename1)
    call mol%print_geometry()
    call mol%read_enuc(filename2)
    print *, "-----------------------------------------------------------------"
    print *, "The overlap matrix is:"
    call mol%read_overlap(filename3)
    ! call dump2(mol%S_C, 7)
    print *, "-----------------------------------------------------------------"
    print *, "The kinetic matrix is:"
    call mol%read_kinetic(filename4)
    ! call dump2(mol%T_C, 7)
    print *, "-----------------------------------------------------------------"
    print *, "The nuclear-attraction matrix is:"
    call mol%read_nuclear_attraction(filename5)
    ! call dump2(mol%V_C, 7)
    print *, "-----------------------------------------------------------------"
    print *, "The core Hamiltonian matrix is:"
    call mol%cons_core_hamiltonian()
    call dump2(mol%H_C, 7)
    print *, "-----------------------------------------------------------------"
    print *, "The two electron integral is:"
    call mol%read_two_electron(filename6)
    ! call dump4(mol%TEI, 7)
    print *, "-----------------------------------------------------------------"
    print *, "The orthogonal matrix is:"
    call mol%cons_orthogonalization_matrix()
    call dump2(mol%S_C_inv, 7)

    ! call mol%cons_core_hamiltonian()
    call mol%cons_initial_guess()
    print *, "-----------------------------------------------------------------"
    print *, "The initial(guess) Fock matrix is:"
    call dump2(mol%F_C, 7)

    call mol%cons_initial_scf_energy()
    print *, "-----------------------------------------------------------------"
    print *, "The initial SCF electronic energy is:", mol%E_elec
    print *, "The total enrgy is:", mol%E_total

    call mol%update_fock_matrix()
    print *, "-----------------------------------------------------------------"
    print *, "The updated Fock matrix is:"
    call dump2(mol%F_C, 7)
    call mol%update_density_matrix()
    print *, "-----------------------------------------------------------------"
    print *, "The updated density matrix is:"
    call dump2(mol%D_C, 7)
    ! call mol%run_scf()
contains

subroutine write_test_file(filename1, filename2, filename3, filename4, filename5, filename6)
    character(len=*), intent(in) :: filename1, filename2, filename3, filename4, filename5, filename6
    integer :: unit

    open(newunit=unit, file=filename1, status='replace', action='write')
    write(unit, '(I1)') 3
    write(unit, '(I1, 3F20.12)') 8,  0.000000000000_dp, -0.143225816552_dp,  0.000000000000_dp
    write(unit, '(I1, 3F20.12)') 1,  1.638036840407_dp,  1.136548822547_dp,  0.000000000000_dp
    write(unit, '(I1, 3F20.12)') 1, -1.638036840407_dp,  1.136548822547_dp,  0.000000000000_dp
    close(unit)

    open(newunit=unit, file=filename2, status='replace', action='write')
    write(unit, '(F20.12)') 8.002367061810450
    close(unit)

    open(newunit=unit, file=filename3, status='replace', action='write')
    write(unit, '(I1, 1X, I1, F20.12)') 1, 1, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1, F20.12)') 2, 1, 0.236703936510848_dp
    write(unit, '(I1, 1X, I1, F20.12)') 2, 2, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1, F20.12)') 3, 1, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1, F20.12)') 3, 2, 0.000000000000000_dp
    write(unit, '(I1, 1X, I1, F20.12)') 3, 3, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 4, 1, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 4, 2, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 4, 3, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 4, 4, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 5, 1, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 5, 2, 0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 5, 3, 0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 5, 4, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 5, 5, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 1, 0.038405599785757_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 2, 0.386138840478249_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 3, 0.268438243716457_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 4, 0.209726941420375_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 5, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 6, 6, 1.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 1, 0.038405599785757_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 2, 0.386138840478250_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 3, -0.268438243716457_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 4, 0.209726941420375_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 5, -0.000000000000000_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 6, 0.181759886298063_dp
    write(unit, '(I1, 1X, I1,F20.12)') 7, 7, 1.000000000000000_dp
    close(unit)

    open(newunit=unit, file=filename4, status='replace', action='write')
    write(unit, '(I1, 1X, I1,F20.12)') 1, 1, 29.003199945539588
    write(unit, '(I1, 1X, I1,F20.12)') 2, 1, -0.168010939316492
    write(unit, '(I1, 1X, I1,F20.12)') 2, 2, 0.808127954930347
    write(unit, '(I1, 1X, I1,F20.12)') 3, 1, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 3, 2, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 3, 3, 2.528731198194763
    write(unit, '(I1, 1X, I1,F20.12)') 4, 1, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 4, 2, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 4, 3, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 4, 4, 2.528731198194763
    write(unit, '(I1, 1X, I1,F20.12)') 5, 1, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 2, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 3, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 4, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 5, 2.528731198194763
    write(unit, '(I1, 1X, I1,F20.12)') 6, 1, -0.008416383544591
    write(unit, '(I1, 1X, I1,F20.12)') 6, 2, 0.070517372751001
    write(unit, '(I1, 1X, I1,F20.12)') 6, 3, 0.147090913304052
    write(unit, '(I1, 1X, I1,F20.12)') 6, 4, 0.114920016354202
    write(unit, '(I1, 1X, I1,F20.12)') 6, 5, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 6, 6, 0.760031883566609
    write(unit, '(I1, 1X, I1,F20.12)') 7, 1, -0.008416383544591
    write(unit, '(I1, 1X, I1,F20.12)') 7, 2, 0.070517372751002
    write(unit, '(I1, 1X, I1,F20.12)') 7, 3, -0.147090913304052
    write(unit, '(I1, 1X, I1,F20.12)') 7, 4, 0.114920016354202
    write(unit, '(I1, 1X, I1,F20.12)') 7, 5, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 7, 6, -0.003979868621841
    write(unit, '(I1, 1X, I1,F20.12)') 7, 7, 0.760031883566609
    close(unit)

    open(newunit=unit, file=filename5, status='replace', action='write')
    write(unit, '(I1, 1X, I1,F20.12)') 1, 1, -61.580595358149914
    write(unit, '(I1, 1X, I1,F20.12)') 2, 1, -7.410821877330996
    write(unit, '(I1, 1X, I1,F20.12)') 2, 2, -10.009071226859687
    write(unit, '(I1, 1X, I1,F20.12)') 3, 1, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 3, 2, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 3, 3, -9.987550507419133
    write(unit, '(I1, 1X, I1,F20.12)') 4, 1, -0.014473835903318
    write(unit, '(I1, 1X, I1,F20.12)') 4, 2, -0.176890246723779
    write(unit, '(I1, 1X, I1,F20.12)') 4, 3, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 4, 4, -9.944042952440439
    write(unit, '(I1, 1X, I1,F20.12)') 5, 1, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 2, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 3, -0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 4, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 5, 5, -9.875876065926128
    write(unit, '(I1, 1X, I1,F20.12)') 6, 1, -1.231685872788711
    write(unit, '(I1, 1X, I1,F20.12)') 6, 2, -2.977227202971528
    write(unit, '(I1, 1X, I1,F20.12)') 6, 3, -1.822241058022723
    write(unit, '(I1, 1X, I1,F20.12)') 6, 4, -1.471788274313766
    write(unit, '(I1, 1X, I1,F20.12)') 6, 5, 0.000000000000000
    write(unit, '(I1, 1X, I1,F20.12)') 6, 6, -5.300202953839792
    write(unit, '(I1, 1X, I1,F20.12)') 7, 1, -1.231685872788712
    write(unit, '(I1, 1X, I1,F20.12)') 7, 2, -2.977227202971529
    write(unit, '(I1, 1X, I1,F20.12)') 7, 3, 1.822241058022724
    write(unit, '(I1, 1X, I1,F20.12)') 7, 4, -1.471788274313766
    write(unit, '(I1, 1X, I1,F20.12)') 7, 5, 0.000000000000001
    write(unit, '(I1, 1X, I1,F20.12)') 7, 6, -1.067166014696110
    write(unit, '(I1, 1X, I1,F20.12)') 7, 7, -5.300202953839793
    close(unit)

    ! open(newunit=unit, file=filename6, status='replace', action='write')
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 1, 1, 1, 1, 4.785065404705_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 2, 1, 1, 1, 0.741380351973_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 2, 2, 1, 1, 1.118946866342_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 2, 1, 2, 1, 0.136873385354_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 2, 2, 2, 1, 0.256633394731_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 2, 2, 2, 2, 0.817206321526_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 3, 1, 1, 1.115813812152_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 4, 1, 1, 1.115813812152_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 1, 1, 1.115813812152_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 1, 3, 1, 0.024477412258_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 1, 4, 1, 0.024477412258_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 1, 5, 1, 0.024477412258_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 3, 2, 1, 0.256683985810_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 4, 2, 1, 0.256683985810_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 2, 1, 0.256683985810_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 2, 3, 1, 0.037808607416_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 2, 4, 1, 0.037808607416_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 2, 5, 1, 0.037808607416_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 3, 2, 2, 0.817022605321_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 4, 2, 2, 0.817022605321_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 2, 2, 0.817022605321_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 2, 3, 2, 0.180518392105_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 2, 4, 2, 0.180518392105_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 2, 5, 2, 0.180518392105_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 3, 3, 3, 3, 0.880159093375_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 3, 4, 3, 0.047444445118_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 4, 3, 3, 0.785270203138_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 4, 4, 4, 4, 0.880159093375_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 3, 5, 3, 0.047444445118_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 4, 5, 4, 0.047444445118_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 3, 3, 0.785270203138_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 4, 4, 0.785270203138_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 5, 5, 5, 5, 0.880159093375_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 1, 1, 0.121785349417_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 1, 1, 0.313334133204_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 2, 1, 0.022309236063_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 2, 1, 0.072840109508_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 2, 2, 0.041611622301_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 2, 2, 0.258884490884_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 1, 1, 0.183538575025_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 1, 1, 0.143396050576_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 3, 1, 0.000807385153_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 4, 1, 0.000630798415_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 2, 1, 0.043197737649_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 2, 1, 0.033749771523_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 3, 1, 0.004317163348_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 4, 1, 0.003372937671_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 3, 2, 0.001519113028_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 4, 2, 0.001186861174_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 2, 2, 0.163910454166_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 2, 2, 0.128060881874_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 3, 2, 0.033774590907_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 2, 4, 2, 0.026387602417_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 3, 1, 0.009644729451_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 4, 1, 0.001744920871_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 3, 1, 0.001744920871_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 4, 1, 0.008774614179_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 5, 5, 1, 0.007411332583_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 3, 3, 0.041685945046_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 4, 3, 0.000121848711_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 4, 4, 0.041625184455_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 1, 5, 5, 0.041529985809_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 3, 2, 0.064012118672_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 3, 4, 2, 0.015709756842_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 3, 2, 0.015709756842_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 4, 4, 2, 0.056178354074_dp
    ! write(unit, '(I1, 1X, I1, 1X, I1, 1X, I1, F20.12)') 6, 5, 
    
    ! close(unit)
end subroutine write_test_file

subroutine dump2(W, N)
    integer, intent(in) :: N
    real(dp), dimension(N, N), intent(in) :: W
    integer :: i
  
    do i = 1, N
      print *, W(i, :)
    end do
end subroutine dump2

subroutine dump4(W, N)
    integer, intent(in) :: N
    real(dp), dimension(N, N, N, N), intent(in) :: W
    integer :: i, j, k, l
  
    do i = 1, N
      do j = 1, N
        do k = 1, N
          do l = 1, N
            print *, i, j, k, l, W(i, j, k, l)
          end do
        end do
      end do
    end do
end subroutine dump4

end program test_mod_scf