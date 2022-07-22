function build_future_df(base_df, to_predict::Symbol)

    subjects = unique(base_df.subject)

    result_lines = Vector{DataFrame}()

    for this_subject in subjects

        subject_df = base_df[ base_df.subject .== this_subject ,:]

        if (size(subject_df, 1) == 1)

            continue;

        else

            for origin_idx in 1:(size(subject_df, 1) - 1)

                if (ismissing(subject_df[origin_idx, :Absiella_dolichum]))

                    continue;

                else

                    for target_idx in (origin_idx + 1):size(subject_df, 1)

                        if (ismissing(subject_df[target_idx, to_predict]))

                            continue;

                        else
                            origin_df = DataFrame()
                            push!(origin_df, copy(subject_df[origin_idx, :]))
                            target_df = subject_df[target_idx, :]

                            insertcols!(origin_df, 1, :target => parse(Float64, target_df[to_predict]))
                            insertcols!(origin_df, 1, :futureAgeMonths => parse(Float64, target_df[:ageMonths]))
                            insertcols!(origin_df, 1, :ageMonthsDelta => parse(Float64, target_df[:ageMonths]) - parse(Float64, origin_df[1, :ageMonths]))
                            insertcols!(origin_df, 1, :futureTimepoint => subject_df[target_idx, :timepoint])
                            insertcols!(origin_df, 1, :timepointDelta => subject_df[target_idx, :timepoint] - subject_df[origin_idx, :timepoint])

                            push!(result_lines, origin_df)

                        end # end if ismissing(to_predict from irigin_idx)

                    end # end for target_idx

                end # end if ismissing(stool from irigin_idx)

            end # end for origin_idx

        end # end if size == 1

    end # end for subject

    return(vcat(result_lines...))

end # end function1