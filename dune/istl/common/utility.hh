// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_ISTL_COMMON_UTILITY_HH
#define DUNE_ISTL_COMMON_UTILITY_HH

#include <optional>

#include <dune/common/parametertree.hh>
#include <dune/common/stdstreams.hh>


namespace Dune::Impl {
  inline int getVerbosity(Dune::ParameterTree param, bool defaultVerbosity=false) {
    std::optional<int> verbosity;
    std::optional<bool> verbose;
    if (param.hasKey("verbose"))
      verbose = param.template get<bool>("verbose");
    if (param.hasKey("verbosity"))
      verbosity = param.template get<int>("verbosity");

    if ((verbose && verbosity) and (verbose.value() != (verbosity.value() > 0)))
      Dune::dwarn <<
        "Both 'verbose' and 'verbosity' are set with conflicting values: "
        "'verbose'=" << (*verbose ? "true" : "false") << ", "
        "'verbosity'=" << *verbosity << ". Please set only one of them." << std::endl;
      if (verbose && !verbosity)
        Dune::dwarn << "Parameter 'verbose' is deprecated. Please use 'verbosity' instead." << std::endl;

    if (verbosity)
      return verbosity.value();
    else
      return verbose.value_or(defaultVerbosity);
  }
} // end namespace Dune::Impl

#endif // DUNE_ISTL_COMMON_UTILITY_HH
