import { IClassificationReview } from "../../models/classification-review.model";
import {
  ILoadReviewAcmgRules,
  ILoadReviewPublications,
} from "../../models/load-review.model";
import { IVUSSummary } from "../../models/vus-summary.model";

export interface ILoadReviewPageRequest {
  vusId: number;
}

export interface ILoadReviewPageResponse {
  variantSummary: IVUSSummary;
  publications: ILoadReviewPublications[];
  acmgRules: ILoadReviewAcmgRules[];
  classifications: string[];
}

export interface ISaveClassificationReviewRequest {
  vusId: number;
  newClassification: string;
  reason?: string;
  publicationIds?: number[];
  acmgRuleIds?: number[];
  isNewAcmgAdded?: boolean;
  isExistingAcmgRemoved?: boolean;
}

export interface ISaveClassificationReviewResponse {
  isSuccess: boolean;
}

export interface IGetAllClassificationReviewsRequest {
  vusId: number;
}

export interface IGetAllClassificationReviewsResponse {
  isSuccess: boolean;
  variantSummary: IVUSSummary;
  reviews: IClassificationReview[];
}
