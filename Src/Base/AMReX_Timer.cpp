// #include <AMReX_Timer.H>

// namespace amrex {

// Timer::Timer ()
// {
//     // m_name = a_name;
//     m_timer_time = 0.0;
//     m_started = false;
//     // if (a_start) {
//     //     start();
//     // }
// }

// void
// Timer::start ()
// {
//     m_start_time = amrex::second();
//     m_started = true;
// }

// void
// Timer::stop ()
// {
//     if (m_started) {
//         m_timer_time += amrex::second() - m_start_time;
//         m_started = false;
//     } else {
//         amrex::Warning("Timer::record: First call Timer::start().\n");
//     }
// }

// void
// Timer::reset ()
// {
//     m_started = false;
//     m_timer_time = 0.0;
// }

// double
// Timer::time ()
// {
//     if (m_started) {
//         stop();
//         start();
//         // add warning about still recording?
//     }
//     return m_timer_time;
// }

// // std::ostream&
// // operator<< (std::ostream& os, Timer& timer)
// // {
// //     os << timer.time();
// //     return os;
// // }


// }