import { useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import ReviewPage from "../components/review-page/review-page";
import { reviewService } from "../services/review/review.service";
import {
  ILoadReviewAcmgRules,
  ILoadReviewPublications,
} from "../models/load-review.model";
import { IVUSSummary } from "../models/vus-summary.model";

const ReviewPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [variantSummary, setVariantSummary] = useState<IVUSSummary>(undefined);
  const [acmgRules, setAcmgRules] = useState<ILoadReviewAcmgRules[]>(undefined);
  const [publications, setPublications] =
    useState<ILoadReviewPublications[]>(undefined);
  const [classifications, setClassifications] = useState<string[]>(undefined);

  const loc = useLocation();
  const vusId = loc.pathname.split("/review/")[1];

  const searchParams = new URLSearchParams(loc.search);
  const acmg = searchParams.get("acmg");
  let selectedAcmgToAdd: ILoadReviewAcmgRules = undefined;

  if (acmg) {
    const acmgSplit = acmg.split("-");
    selectedAcmgToAdd = { id: parseInt(acmgSplit[0]), name: acmgSplit[1] };
  }

  const acmgRemove = searchParams.get("acmg-remove");
  let selectedAcmgToRemove: ILoadReviewAcmgRules = undefined;

  if (acmgRemove) {
    const acmgSplit = acmgRemove.split("-");
    selectedAcmgToRemove = { id: parseInt(acmgSplit[0]), name: acmgSplit[1] };
  }

  useEffect(() => {
    if (isLoading) {
      reviewService.loadReviewPage({ vusId: parseInt(vusId) }).then((res) => {
        setVariantSummary(res.variantSummary);
        setAcmgRules(res.acmgRules);
        setPublications(res.publications);
        setClassifications(res.classifications);
        setIsLoading(false);
      });
    }
  }, [isLoading, vusId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <ReviewPage
        variantSummary={variantSummary}
        acmgRules={acmgRules}
        publications={publications}
        classifications={classifications}
        selectedAcmgRuleToAdd={selectedAcmgToAdd}
        selectedAcmgRuleToRemove={selectedAcmgToRemove}
        reviewService={reviewService}
      />
    );
  }
};

export default ReviewPageWrapper;
